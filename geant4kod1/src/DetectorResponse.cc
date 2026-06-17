#include "DetectorResponse.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>

DetectorResponse::DetectorResponse()
    : fElectricField(0.273),
      fDriftVelocity(0.1098),
      fElectronLifetime(18000.0),
      fSamplingPeriod(0.5),
      fWireSpacing(0.3),
      fLArDensity(1.396),
      fTemperature(89.0),
      fWionization(23.6e-6),

      fRecombAlpha(0.930),
      fRecombBeta(0.212),

      fDiffusionL(6.2),

      fDiffusionT(9.1),

      fENC(400.0),
      fADCPerElectron(1.0/6.8),
      fShapingTime(2.0),

      fFieldResponseTime(3.0),
      fResponseNTicks(20),

      fNoise1fKnee(10000.0),
      fCoherentNoiseRMS(100.0),

      fWienerFilterCutoff(150000.0)
{
}

DetectorResponse::~DetectorResponse()
{}

G4double DetectorResponse::CalculateRecombination(G4double dEdx_MeV_per_cm) const
{
    if (dEdx_MeV_per_cm <= 0.0) return 0.0;

    G4double beta_eff = fRecombBeta / (fLArDensity * fElectricField);

    G4double betadEdx = beta_eff * dEdx_MeV_per_cm;

    G4double R = std::log(fRecombAlpha + betadEdx) / betadEdx;

    R = std::max(0.0, std::min(1.0, R));

    return R;
}

G4double DetectorResponse::EnergyToElectrons(G4double energyMeV,
                                              G4double dEdx_MeV_per_cm) const
{
    if (energyMeV <= 0.0) return 0.0;

    G4double n_total = energyMeV / fWionization;

    G4double R = CalculateRecombination(dEdx_MeV_per_cm);
    G4double n_surviving = n_total * R;

    return n_surviving;
}

G4double DetectorResponse::CalculateDriftTime(G4double distanceToPlaneCm) const
{
    if (distanceToPlaneCm <= 0.0) return 0.0;
    return distanceToPlaneCm / fDriftVelocity;
}

G4double DetectorResponse::CalculateLifetimeAttenuation(G4double driftTimeMicrosec) const
{
    if (driftTimeMicrosec <= 0.0) return 1.0;
    return std::exp(-driftTimeMicrosec / fElectronLifetime);
}

G4double DetectorResponse::CalculateLongitudinalDiffusion(G4double driftTimeMicrosec) const
{
    if (driftTimeMicrosec <= 0.0) return 0.0;

    G4double t_sec = driftTimeMicrosec * 1e-6;

    G4double sigma_cm = std::sqrt(2.0 * fDiffusionL * t_sec);

    G4double sigma_ticks = sigma_cm / (fDriftVelocity * fSamplingPeriod);

    return sigma_ticks;
}

G4double DetectorResponse::CalculateTransverseDiffusion(G4double driftTimeMicrosec) const
{
    if (driftTimeMicrosec <= 0.0) return 0.0;

    G4double t_sec = driftTimeMicrosec * 1e-6;

    G4double sigma_cm = std::sqrt(2.0 * fDiffusionT * t_sec);

    G4double sigma_wires = sigma_cm / fWireSpacing;

    return sigma_wires;
}

G4double DetectorResponse::GetElectronicsNoise() const
{
    return G4RandGauss::shoot(0.0, fENC);
}

G4double DetectorResponse::GetADCConversion() const
{
    return fADCPerElectron;
}

std::vector<G4double> DetectorResponse::GetFieldResponse(G4int planeID) const
{
    std::vector<G4double> response(fResponseNTicks, 0.0);

    for (G4int i = 0; i < fResponseNTicks; i++) {
        G4double t = i * fSamplingPeriod;
        G4double tNorm = t / fFieldResponseTime;

        if (planeID == 2) {
            response[i] = tNorm * std::exp(1.0 - tNorm);
        } else {
            response[i] = (1.0 - tNorm) * std::exp(1.0 - tNorm);
        }
    }

    G4double norm = 0.0;
    if (planeID == 2) {
        for (auto v : response) norm += v;
        if (norm > 0) for (auto& v : response) v /= norm;
    } else {
        for (auto v : response) norm = std::max(norm, std::abs(v));
        if (norm > 0) for (auto& v : response) v /= norm;
    }

    return response;
}

std::vector<G4double> DetectorResponse::GetElectronicsResponse() const
{
    std::vector<G4double> response(fResponseNTicks, 0.0);

    G4double norm = 0.0;
    for (G4int i = 0; i < fResponseNTicks; i++) {
        G4double t = i * fSamplingPeriod;
        G4double tNorm = t / fShapingTime;

        G4double tN4 = tNorm * tNorm * tNorm * tNorm;
        response[i] = tN4 * std::exp(-4.0 * tNorm);
        norm += response[i];
    }

    if (norm > 0.0) {
        for (auto& v : response) v /= norm;
    }

    return response;
}

std::vector<G4double> DetectorResponse::GetSignalResponse(G4int planeID) const
{
    std::vector<G4double> fieldResp = GetFieldResponse(planeID);
    std::vector<G4double> elecResp = GetElectronicsResponse();

    G4int nField = fieldResp.size();
    G4int nElec = elecResp.size();
    G4int nConv = nField + nElec - 1;

    std::vector<G4double> combined(nConv, 0.0);

    for (G4int i = 0; i < nConv; i++) {
        G4double sum = 0.0;
        for (G4int j = 0; j < nField; j++) {
            G4int k = i - j;
            if (k >= 0 && k < nElec) {
                sum += fieldResp[j] * elecResp[k];
            }
        }
        combined[i] = sum;
    }

    return combined;
}

const DetectorResponse::TrigTable& DetectorResponse::GetCachedTrigTable(G4int N)
{
    auto it = fTrigTableCache.find(N);
    if (it != fTrigTableCache.end()) return it->second;

    TrigTable& table = fTrigTableCache[N];
    table.cosTable.resize(N);
    table.sinTable.resize(N);
    for (G4int k = 0; k < N; k++) {
        G4double angle = 2.0 * M_PI * k / N;
        table.cosTable[k] = std::cos(angle);
        table.sinTable[k] = std::sin(angle);
    }
    return table;
}

const DetectorResponse::ResponseDFTCache& DetectorResponse::GetCachedResponseDFT(
    const std::vector<G4double>& signalResponse, G4int planeID, G4int N)
{
    G4int key = planeID * 100000 + N;
    auto it = fResponseDFTCache.find(key);
    if (it != fResponseDFTCache.end()) return it->second;

    G4int nFreq = N / 2 + 1;
    G4int nResp = signalResponse.size();

    ResponseDFTCache& cache = fResponseDFTCache[key];
    cache.re.resize(nFreq, 0.0);
    cache.im.resize(nFreq, 0.0);

    const TrigTable& trig = GetCachedTrigTable(N);

    for (G4int k = 0; k < nFreq; k++) {
        G4double sumRe = 0.0, sumIm = 0.0;
        for (G4int n = 0; n < nResp && n < N; n++) {
            G4int idx = (static_cast<G4long>(k) * n) % N;
            sumRe += signalResponse[n] * trig.cosTable[idx];
            sumIm -= signalResponse[n] * trig.sinTable[idx];
        }
        cache.re[k] = sumRe;
        cache.im[k] = sumIm;
    }

    return cache;
}

void DetectorResponse::GenerateNoise(std::vector<float>& waveform, G4int planeID) const
{
    G4int N = waveform.size();
    if (N < 2) return;

    G4double samplingRate = 1.0 / fSamplingPeriod;
    G4double samplingRateHz = samplingRate * 1e6;
    G4double df = samplingRateHz / N;

    G4int nFreq = N / 2 + 1;

    std::vector<G4double> ampSpectrum(nFreq, 0.0);
    std::vector<G4double> phases(nFreq, 0.0);

    G4double whiteAmp = fENC * std::sqrt(2.0 / N);

    for (G4int k = 0; k < nFreq; k++) {
        G4double f = k * df;
        if (f < 1.0) f = 1.0;

        G4double intrinsicPower = 1.0 + fNoise1fKnee / f;

        G4double coherentPower = 0.0;
        if (f < 50000.0) {
            G4double coherentAmp = fCoherentNoiseRMS * std::sqrt(2.0 / N);
            coherentPower = coherentAmp * coherentAmp * std::exp(-f / 20000.0);
        }

        G4double totalAmp = whiteAmp * std::sqrt(intrinsicPower) +
                            std::sqrt(coherentPower);

        if (planeID == 0) totalAmp *= 1.15;
        else if (planeID == 1) totalAmp *= 1.08;

        ampSpectrum[k] = totalAmp;
        phases[k] = G4UniformRand() * 2.0 * M_PI;
    }

    ampSpectrum[0] = 0.0;

    std::vector<G4double> specRe(nFreq), specIm(nFreq);
    for (G4int k = 0; k < nFreq; k++) {
        specRe[k] = ampSpectrum[k] * std::cos(phases[k]);
        specIm[k] = ampSpectrum[k] * std::sin(phases[k]);
    }

    G4double scale = 2.0 / std::sqrt(static_cast<double>(N));

    std::vector<G4double> cosBase(N), sinBase(N);
    for (G4int i = 0; i < N; i++) {
        G4double angle = 2.0 * M_PI * i / N;
        cosBase[i] = std::cos(angle);
        sinBase[i] = std::sin(angle);
    }

    for (G4int n = 0; n < N; n++) {
        G4double val = 0.0;
        for (G4int k = 1; k < nFreq - 1; k++) {
            G4int idx = (static_cast<G4long>(k) * n) % N;
            val += specRe[k] * cosBase[idx] - specIm[k] * sinBase[idx];
        }
        if (N % 2 == 0 && nFreq > 1) {
            G4double cosNyq = (n % 2 == 0) ? 1.0 : -1.0;
            val += 0.5 * ampSpectrum[nFreq-1] * std::cos(phases[nFreq-1]) * cosNyq;
        }
        waveform[n] += static_cast<float>(val * scale);
    }
}

std::vector<float> DetectorResponse::ApplyWienerDeconvolution(
    const std::vector<float>& rawWaveform,
    const std::vector<G4double>& signalResponse,
    G4double noiseRMS,
    G4int planeID)
{
    G4int N = rawWaveform.size();
    std::vector<float> result(N, 0.0f);
    if (N < 4) return result;

    G4int nFreq = N / 2 + 1;
    G4double samplingRateHz = 1e6 / fSamplingPeriod;
    G4double df = samplingRateHz / N;

    const TrigTable& trig = GetCachedTrigTable(N);

    std::vector<G4double> rawRe(nFreq, 0.0), rawIm(nFreq, 0.0);
    for (G4int k = 0; k < nFreq; k++) {
        G4double sumRe = 0.0, sumIm = 0.0;
        for (G4int n = 0; n < N; n++) {
            G4int idx = (static_cast<G4long>(k) * n) % N;
            sumRe += rawWaveform[n] * trig.cosTable[idx];
            sumIm -= rawWaveform[n] * trig.sinTable[idx];
        }
        rawRe[k] = sumRe;
        rawIm[k] = sumIm;
    }

    const ResponseDFTCache& respDFT = GetCachedResponseDFT(signalResponse, planeID, N);
    const std::vector<G4double>& respRe = respDFT.re;
    const std::vector<G4double>& respIm = respDFT.im;

    G4double signalPower = 0.0;
    for (G4int k = 0; k < nFreq; k++) {
        signalPower += rawRe[k]*rawRe[k] + rawIm[k]*rawIm[k];
    }
    signalPower /= nFreq;
    G4double noisePower = noiseRMS * noiseRMS * N;
    G4double NSR = (signalPower > 0) ? noisePower / signalPower : 0.01;
    NSR = std::max(0.001, std::min(1.0, NSR));

    std::vector<G4double> deconvRe(nFreq, 0.0), deconvIm(nFreq, 0.0);

    for (G4int k = 0; k < nFreq; k++) {
        G4double f = k * df;

        G4double Hpower = respRe[k]*respRe[k] + respIm[k]*respIm[k];

        G4double denom = Hpower + NSR;
        if (denom < 1e-20) denom = 1e-20;

        G4double filterRe = respRe[k] / denom;
        G4double filterIm = -respIm[k] / denom;

        G4double sRe = rawRe[k]*filterRe - rawIm[k]*filterIm;
        G4double sIm = rawRe[k]*filterIm + rawIm[k]*filterRe;

        G4double gaussianFilter = std::exp(-f*f / (2.0 * fWienerFilterCutoff * fWienerFilterCutoff));

        deconvRe[k] = sRe * gaussianFilter;
        deconvIm[k] = sIm * gaussianFilter;
    }

    for (G4int n = 0; n < N; n++) {
        G4double val = deconvRe[0];
        for (G4int k = 1; k < nFreq - 1; k++) {
            G4int idx = (static_cast<G4long>(k) * n) % N;
            val += 2.0 * (deconvRe[k] * trig.cosTable[idx] - deconvIm[k] * trig.sinTable[idx]);
        }
        if (N % 2 == 0 && nFreq > 1) {
            G4double cosNyq = (n % 2 == 0) ? 1.0 : -1.0;
            val += deconvRe[nFreq-1] * cosNyq;
        }
        result[n] = static_cast<float>(val / N);
    }

    return result;
}

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "EventInformation.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionZero.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4TauMinus.hh"

#include <cmath>
#include <vector>

namespace {
    constexpr G4double kMicroBooNEReadoutWindow = 4.8 * millisecond;
    constexpr G4double kMicroBooNEBeamSpillTime = 1.6 * millisecond;

    constexpr G4double kDefaultCosmicOverlayMean = 8.0;

    G4double GetMaxVisiblePrimaryTime()
    {
        const G4double maxDriftTimeUs =
            DetectorConstruction::fTPC_X / DetectorConstruction::fDriftVelocity;
        const G4double maxVisible = kMicroBooNEReadoutWindow - maxDriftTimeUs * microsecond;
        return (maxVisible > 0.0) ? maxVisible : 0.0;
    }
}

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(nullptr),
      fGeneratorMode(kNeutrinoMode),
      fNeutrinoFlavor(kNuMu),
      fSamplingProfile(kRealisticProfile),
      fCosmicMinEnergy(0.5*GeV),
      fCosmicMaxEnergy(1000.*GeV),
      fCosmicMinTheta(0.*deg),
      fCosmicMaxTheta(60.*deg),
      fNeutrinoMinEnergy(0.2*GeV),
      fNeutrinoMaxEnergy(5.*GeV),
      fCosmicEventCount(0),
      fNeutrinoEventCount(0),
      fNuMuCCCount(0),
      fNueCCCount(0),
      fNuTauCCCount(0),
      fNCCount(0),
      f2p2hCount(0),
      fFSIAbsorptionCount(0),
      fFSIChargeExchangeCount(0),
      fNuMuFraction(0.25),
      fNuEFraction(0.25),
      fNuTauFraction(0.25),
      fNCFraction(0.25),
      fEnableCosmicOverlay(true),
      fCosmicOverlayMean(kDefaultCosmicOverlayMean)
{
    fParticleGun = new G4ParticleGun(1);

    fParticleGun->SetParticleDefinition(G4MuonMinus::MuonMinus());
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
    fParticleGun->SetParticleEnergy(4.*GeV);

    InitializeBNBFlux();

    fMessenger = new PrimaryGeneratorMessenger(this);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
    delete fMessenger;

    G4cout << "Event generation summary: "
           << "cosmic=" << fCosmicEventCount
           << " neutrino=" << fNeutrinoEventCount
           << " (numu=" << fNuMuCCCount
           << " nue=" << fNueCCCount
           << " nutau=" << fNuTauCCCount
           << " nc=" << fNCCount
           << " 2p2h=" << f2p2hCount << ")"
           << " FSI[abs=" << fFSIAbsorptionCount
           << " cex=" << fFSIChargeExchangeCount << "]" << G4endl;
}

void PrimaryGeneratorAction::SetFlavorRatios(G4double nuMuFrac, G4double nuEFrac,
                                              G4double nuTauFrac, G4double ncFrac)
{
    G4double total = nuMuFrac + nuEFrac + nuTauFrac + ncFrac;
    if (total <= 0.0) return;
    fNuMuFraction = nuMuFrac / total;
    fNuEFraction = nuEFrac / total;
    fNuTauFraction = nuTauFrac / total;
    fNCFraction = ncFrac / total;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    EventInformation* eventInfo = new EventInformation();
    anEvent->SetUserInformation(eventInfo);

    switch(fGeneratorMode) {
        case kCosmicRayMode:
            GenerateCosmicRay(anEvent, eventInfo);
            fCosmicEventCount++;
            break;
        case kNeutrinoMode:
            GenerateNeutrinoEvent(anEvent, eventInfo);
            fNeutrinoEventCount++;
            break;
        case kTestMode:
            GenerateTestParticle(anEvent);
            break;
    }
}

void PrimaryGeneratorAction::GenerateNeutrinoEvent(G4Event* anEvent, EventInformation* eventInfo)
{
    G4double Enu = SampleBNBEnergy();
    fParticleGun->SetParticleTime(kMicroBooNEBeamSpillTime);

    NeutrinoFlavor flavor = fNeutrinoFlavor;

    if (flavor == kAllFlavors) {
        if (fSamplingProfile == kMLBalancedProfile) {
            G4double rand = G4UniformRand();
            if (rand < fNuMuFraction) {
                flavor = kNuMu;
            } else if (rand < fNuMuFraction + fNuEFraction) {
                flavor = kNuE;
            } else if (rand < fNuMuFraction + fNuEFraction + fNuTauFraction) {
                flavor = kNuTau;
            } else {
                flavor = kNC;
            }
        } else {
            flavor = kNuMu;
        }
    }

    switch (flavor) {
        case kNuMu: {
            G4double rand = G4UniformRand();
            if (rand < 0.35) {
                GenerateNuMuCCQE(anEvent, eventInfo, Enu);
            } else if (rand < 0.50) {
                Generate2p2hMEC(anEvent, eventInfo, Enu, 13);
                f2p2hCount++;
            } else if (rand < 0.72) {
                GenerateNuMuCC1Pi(anEvent, eventInfo, Enu);
            } else {
                GenerateNuMuDIS(anEvent, eventInfo, Enu);
            }
            fNuMuCCCount++;
            break;
        }
        case kNuE: {
            G4double rand = G4UniformRand();
            if (rand < 0.35) {
                GenerateNueCCQE(anEvent, eventInfo, Enu);
            } else if (rand < 0.50) {
                Generate2p2hMEC(anEvent, eventInfo, Enu, 11);
                f2p2hCount++;
            } else if (rand < 0.72) {
                GenerateNueCC1Pi(anEvent, eventInfo, Enu);
            } else {
                GenerateNueDIS(anEvent, eventInfo, Enu);
            }
            fNueCCCount++;
            break;
        }
        case kNuTau: {
            if (fSamplingProfile == kRealisticProfile) {
                GenerateNuMuCCQE(anEvent, eventInfo, Enu);
                fNuMuCCCount++;
                break;
            }
            G4double tauThreshold = 3.5 * GeV;
            if (Enu < tauThreshold) {
                Enu = tauThreshold + G4UniformRand() * 6.5 * GeV;
            }

            G4double rand = G4UniformRand();
            if (rand < 0.50) {
                GenerateNuTauCCQE(anEvent, eventInfo, Enu);
            } else {
                GenerateNuTauDIS(anEvent, eventInfo, Enu);
            }
            fNuTauCCCount++;
            break;
        }
        case kNC: {
            G4double rand = G4UniformRand();
            if (rand < 0.40) {
                GenerateNCQE(anEvent, eventInfo, Enu);
            } else if (rand < 0.70) {
                GenerateNCRes(anEvent, eventInfo, Enu);
            } else {
                GenerateNCDIS(anEvent, eventInfo, Enu);
            }
            fNCCount++;
            break;
        }
        default:
            GenerateNuMuCCQE(anEvent, eventInfo, Enu);
            fNuMuCCCount++;
            break;
    }

    if (fSamplingProfile == kRealisticProfile && fEnableCosmicOverlay) {
        GenerateCosmicOverlay(anEvent);
    }
}

G4ThreeVector PrimaryGeneratorAction::SampleVertexInTPC()
{
    G4double x = (G4UniformRand() - 0.5) * DetectorConstruction::fTPC_X * cm;
    G4double y = (G4UniformRand() - 0.5) * DetectorConstruction::fTPC_Y * cm;
    G4double z = (G4UniformRand() - 0.5) * DetectorConstruction::fTPC_Z * cm;
    return G4ThreeVector(x, y, z);
}

void PrimaryGeneratorAction::GenerateNuMuCCQE(G4Event* anEvent, EventInformation* eventInfo, G4double Enu)
{
    G4ThreeVector position = SampleVertexInTPC();

    G4double m_mu = 0.1057*GeV;
    G4double M_n = 0.939*GeV;
    G4double M_p = 0.938*GeV;

    G4double Q2 = SampleCCQE_Q2(Enu);

    G4double kF = 0.220*GeV;
    G4double p_fermi = kF * std::cbrt(G4UniformRand());
    G4ThreeVector pF_vec = G4RandomDirection() * p_fermi;
    G4double EF = sqrt(M_n*M_n + p_fermi*p_fermi);
    G4double Eb = 27*MeV;

    G4double Eavail = Enu + EF - Eb;
    G4double Emu = Eavail - Q2/(2.0*M_n);

    if (Emu < m_mu + 50*MeV) Emu = m_mu + 50*MeV;
    if (Emu > Enu - 50*MeV) Emu = Enu - 50*MeV;

    G4double cos_theta = 1.0 - (Q2 + m_mu*m_mu)/(2.0*Enu*Emu);
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

    G4double theta_mu = acos(cos_theta);
    G4double phi_mu = G4UniformRand() * 2.0 * pi;

    G4ThreeVector mu_dir(sin(theta_mu)*cos(phi_mu),
                         sin(theta_mu)*sin(phi_mu),
                         cos(theta_mu));

    G4double KEmu = Emu - m_mu;
    if (KEmu <= 0) KEmu = 50*MeV;

    fParticleGun->SetParticleDefinition(G4MuonMinus::MuonMinus());
    fParticleGun->SetParticleEnergy(KEmu);
    fParticleGun->SetParticleMomentumDirection(mu_dir);
    fParticleGun->SetParticlePosition(position);
    fParticleGun->GeneratePrimaryVertex(anEvent);

    G4double p_mu = (Emu*Emu - m_mu*m_mu > 0) ? sqrt(Emu*Emu - m_mu*m_mu) : 0.05*GeV;
    G4ThreeVector p_mu_vec = mu_dir * p_mu;
    G4ThreeVector p_nu_vec(0, 0, Enu);
    G4ThreeVector p_p_vec = p_nu_vec + pF_vec - p_mu_vec;

    G4double p_p = p_p_vec.mag();
    G4double Ep = sqrt(M_p*M_p + p_p*p_p);
    G4double KEp = Ep - M_p;

    if (KEp > 10*MeV && p_p > 0 && PassesPauliBlocking(p_p)) {
        fParticleGun->SetParticleDefinition(G4Proton::Proton());
        fParticleGun->SetParticleEnergy(KEp);
        fParticleGun->SetParticleMomentumDirection(p_p_vec.unit());
        fParticleGun->SetParticlePosition(position);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    eventInfo->SetPrimaryVertex(position);
    eventInfo->SetPrimaryMomentum(p_mu_vec/GeV);
    eventInfo->SetPrimaryEnergy(Emu / GeV);
    eventInfo->SetPrimaryPDG(13);

    G4double W2_val = M_n*M_n + 2*M_n*(Enu-Emu) - Q2;
    G4double W = (W2_val > 0.01*GeV*GeV) ? sqrt(W2_val) : 0.939*GeV;

    eventInfo->SetNeutrinoKinematics(Q2/(GeV*GeV), W/GeV, Enu/GeV, "NuMu_CCQE");
}

void PrimaryGeneratorAction::GenerateNuMuCC1Pi(G4Event* anEvent, EventInformation* eventInfo, G4double Enu)
{
    if (Enu < 0.35*GeV) {
        GenerateNuMuCCQE(anEvent, eventInfo, Enu);
        return;
    }

    G4ThreeVector position = SampleVertexInTPC();

    G4double m_mu = 0.1057*GeV;
    G4double M_N = 0.939*GeV;
    G4double W = SampleResonance_W();

    G4double Q2_scale = 0.5*GeV*GeV;
    G4double Q2 = -Q2_scale * log(G4UniformRand());
    G4double Q2_max = 2.0*M_N*Enu;
    if (Q2 > Q2_max) Q2 = Q2_max * G4UniformRand();
    if (Q2 < 0.02*GeV*GeV) Q2 = 0.02*GeV*GeV;

    G4double nu = (W*W - M_N*M_N + Q2) / (2.0*M_N);
    G4double Emu = Enu - nu;

    if (Emu < m_mu + 80*MeV) Emu = m_mu + 80*MeV;
    if (Emu > Enu - 200*MeV) Emu = Enu - 200*MeV;

    G4double cos_theta = 1.0 - (Q2 + m_mu*m_mu)/(2.0*Enu*Emu);
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

    G4double theta_mu = acos(cos_theta);
    G4double phi_mu = G4UniformRand() * 2.0 * pi;

    G4ThreeVector mu_dir(sin(theta_mu)*cos(phi_mu),
                         sin(theta_mu)*sin(phi_mu),
                         cos(theta_mu));

    G4double KEmu = Emu - m_mu;
    if (KEmu <= 0) KEmu = 80*MeV;

    fParticleGun->SetParticleDefinition(G4MuonMinus::MuonMinus());
    fParticleGun->SetParticleEnergy(KEmu);
    fParticleGun->SetParticleMomentumDirection(mu_dir);
    fParticleGun->SetParticlePosition(position);
    fParticleGun->GeneratePrimaryVertex(anEvent);

    G4double p_mu = (Emu*Emu - m_mu*m_mu > 0) ? sqrt(Emu*Emu - m_mu*m_mu) : 0.05*GeV;
    G4ThreeVector p_mu_vec = mu_dir * p_mu;

    G4double E_available = Enu - Emu;
    G4double E_remaining = E_available;

    G4double min_p = 80*MeV, min_pi = 80*MeV;
    if (E_remaining < min_p + min_pi) {
        G4double scale = E_remaining / (min_p + min_pi);
        min_p *= scale * 0.6;
        min_pi *= scale * 0.4;
    }

    G4double KEp = min_p + (E_remaining - min_p - min_pi) * (0.4 + 0.3*G4UniformRand());
    if (KEp > E_remaining - min_pi - 10*MeV) KEp = E_remaining - min_pi - 10*MeV;
    E_remaining -= KEp;

    G4double KEpi = E_remaining * (0.6 + 0.4*G4UniformRand());
    if (KEpi > E_remaining - 10*MeV) KEpi = E_remaining - 10*MeV;

    G4int pionPDG;
    G4double randPion = G4UniformRand();
    if (randPion < 0.55) pionPDG = 211;
    else                 pionPDG = 111;

    G4ThreeVector p_nu_vec_cc1pi(0, 0, Enu);
    G4ThreeVector q_vec_cc1pi = p_nu_vec_cc1pi - p_mu_vec;
    G4ThreeVector q_dir_cc1pi = (q_vec_cc1pi.mag() > 0) ? q_vec_cc1pi.unit()
                                                        : G4ThreeVector(0,0,1);

    std::vector<G4int> hadPDGs;
    std::vector<G4double> hadKEs;
    std::vector<G4ThreeVector> hadDirs;

    if (KEp > 20*MeV) {
        hadPDGs.push_back(2212);
        hadKEs.push_back(KEp);
        hadDirs.push_back(SampleForwardDirection(q_dir_cc1pi, 0.35));
    }
    if (KEpi > 30*MeV) {
        hadPDGs.push_back(pionPDG);
        hadKEs.push_back(KEpi);
        hadDirs.push_back(SampleForwardDirection(q_dir_cc1pi, 0.45));
    }

    ApplyFSI(anEvent, position, hadPDGs, hadKEs, hadDirs);

    eventInfo->SetPrimaryVertex(position);
    eventInfo->SetPrimaryMomentum(p_mu_vec/GeV);
    eventInfo->SetPrimaryEnergy(Emu / GeV);
    eventInfo->SetPrimaryPDG(13);
    eventInfo->SetNeutrinoKinematics(Q2/(GeV*GeV), W/GeV, Enu/GeV, "NuMu_CC1Pi");
}

void PrimaryGeneratorAction::GenerateNuMuDIS(G4Event* anEvent, EventInformation* eventInfo, G4double Enu)
{
    if (Enu < 0.5*GeV) {
        GenerateNuMuCCQE(anEvent, eventInfo, Enu);
        return;
    }

    G4ThreeVector position = SampleVertexInTPC();

    G4double m_mu = 0.1057*GeV;
    G4double M_N = 0.939*GeV;

    G4double y = SampleDIS_y(Enu);
    G4double nu = y * Enu;
    G4double Emu = Enu - nu;

    if (Emu < m_mu + 100*MeV) {
        Emu = m_mu + 100*MeV;
        nu = Enu - Emu;
    }

    G4double Q2;
    if (G4UniformRand() < 0.7)
        Q2 = -2.0*GeV*GeV * log(G4UniformRand());
    else
        Q2 = 0.5*GeV*GeV + 2.5*GeV*GeV * G4UniformRand();

    G4double Q2_max = 2.0*M_N*Enu;
    if (Q2 > Q2_max) Q2 = Q2_max * 0.8;
    if (Q2 < 0.1*GeV*GeV) Q2 = 0.1*GeV*GeV;

    G4double W2 = M_N*M_N + 2.0*M_N*nu - Q2;
    G4double W;
    if (W2 > 0) {
        W = sqrt(W2);
        if (W < 1.1*GeV) W = 1.1*GeV + 0.5*GeV * G4UniformRand();
    } else {
        W = 1.1*GeV + 0.3*GeV * G4UniformRand();
    }

    G4double cos_theta = 1.0 - (Q2 + m_mu*m_mu)/(2.0*Enu*Emu);
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

    G4double theta_mu = acos(cos_theta);
    G4double phi_mu = G4UniformRand() * 2.0 * pi;

    G4ThreeVector mu_dir(sin(theta_mu)*cos(phi_mu),
                         sin(theta_mu)*sin(phi_mu),
                         cos(theta_mu));

    G4double KEmu = Emu - m_mu;
    if (KEmu <= 0) KEmu = 100*MeV;

    fParticleGun->SetParticleDefinition(G4MuonMinus::MuonMinus());
    fParticleGun->SetParticleEnergy(KEmu);
    fParticleGun->SetParticleMomentumDirection(mu_dir);
    fParticleGun->SetParticlePosition(position);
    fParticleGun->GeneratePrimaryVertex(anEvent);

    G4double p_mu = (Emu*Emu - m_mu*m_mu > 0) ? sqrt(Emu*Emu - m_mu*m_mu) : 0.05*GeV;
    G4ThreeVector p_mu_vec = mu_dir * p_mu;

    G4double E_available = Enu - Emu;
    G4int nParticles = 5 + G4int(G4UniformRand() * 5);

    std::vector<G4ParticleDefinition*> hadrons = {
        G4Proton::Proton(), G4Neutron::Neutron(),
        G4PionPlus::PionPlus(), G4PionMinus::PionMinus(),
        G4PionZero::PionZero()
    };

    G4ThreeVector p_nu_vec_dis(0, 0, Enu);
    G4ThreeVector q_vec = p_nu_vec_dis - p_mu_vec;
    G4ThreeVector q_dir = (q_vec.mag() > 0) ? q_vec.unit() : G4ThreeVector(0,0,1);

    G4double E_per = E_available / nParticles;
    G4double E_rem = E_available;

    for (G4int i = 0; i < nParticles; i++) {
        G4double frac = exp(-2.5 * G4UniformRand());
        G4double E_had = E_per * frac * 2.5;

        if (E_had < 70*MeV) E_had = 70*MeV;
        if (E_had > E_rem * 0.7) E_had = E_rem * 0.7;
        if (i == nParticles - 1) E_had = E_rem * 0.8;

        E_rem -= E_had;
        if (E_rem < 0) E_rem = 0;

        if (E_had > 60*MeV) {

            G4ParticleDefinition* hadron;
            G4double rH = G4UniformRand();
            if      (rH < 0.22) hadron = G4Proton::Proton();
            else if (rH < 0.44) hadron = G4Neutron::Neutron();
            else if (rH < 0.63) hadron = G4PionPlus::PionPlus();
            else if (rH < 0.82) hadron = G4PionMinus::PionMinus();
            else                hadron = G4PionZero::PionZero();

            G4double thetaScale = (i == 0) ? 0.25 : 0.40;
            G4ThreeVector had_dir = SampleForwardDirection(q_dir, thetaScale);

            fParticleGun->SetParticleDefinition(hadron);
            fParticleGun->SetParticleEnergy(E_had);
            fParticleGun->SetParticleMomentumDirection(had_dir);
            fParticleGun->SetParticlePosition(position);
            fParticleGun->GeneratePrimaryVertex(anEvent);
        }
    }

    eventInfo->SetPrimaryVertex(position);
    eventInfo->SetPrimaryMomentum(p_mu_vec/GeV);
    eventInfo->SetPrimaryEnergy(Emu / GeV);
    eventInfo->SetPrimaryPDG(13);
    eventInfo->SetNeutrinoKinematics(Q2/(GeV*GeV), W/GeV, Enu/GeV, "NuMu_DIS");
}

void PrimaryGeneratorAction::GenerateNueCCQE(G4Event* anEvent, EventInformation* eventInfo, G4double Enu)
{
    G4ThreeVector position = SampleVertexInTPC();

    G4double m_e = 0.000511*GeV;
    G4double M_n = 0.939*GeV;
    G4double M_p = 0.938*GeV;

    G4double Q2 = SampleCCQE_Q2(Enu);

    G4double kF = 0.220*GeV;
    G4double p_fermi = kF * std::cbrt(G4UniformRand());
    G4ThreeVector pF_vec = G4RandomDirection() * p_fermi;
    G4double EF = sqrt(M_n*M_n + p_fermi*p_fermi);
    G4double Eb = 27*MeV;

    G4double Eavail = Enu + EF - Eb;
    G4double Ee = Eavail - Q2/(2.0*M_n);

    if (Ee < m_e + 10*MeV) Ee = m_e + 10*MeV;
    if (Ee > Enu - 30*MeV) Ee = Enu - 30*MeV;

    G4double cos_theta = 1.0 - (Q2 + m_e*m_e)/(2.0*Enu*Ee);
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

    G4double theta_e = acos(cos_theta);
    G4double phi_e = G4UniformRand() * 2.0 * pi;

    G4ThreeVector e_dir(sin(theta_e)*cos(phi_e),
                        sin(theta_e)*sin(phi_e),
                        cos(theta_e));

    G4double KEe = Ee - m_e;
    if (KEe <= 0) KEe = 10*MeV;

    fParticleGun->SetParticleDefinition(G4Electron::Electron());
    fParticleGun->SetParticleEnergy(KEe);
    fParticleGun->SetParticleMomentumDirection(e_dir);
    fParticleGun->SetParticlePosition(position);
    fParticleGun->GeneratePrimaryVertex(anEvent);

    G4double p_e = (Ee*Ee - m_e*m_e > 0) ? sqrt(Ee*Ee - m_e*m_e) : 0.01*GeV;
    G4ThreeVector p_e_vec = e_dir * p_e;
    G4ThreeVector p_nu_vec(0, 0, Enu);
    G4ThreeVector p_p_vec = p_nu_vec + pF_vec - p_e_vec;

    G4double p_p = p_p_vec.mag();
    G4double Ep = sqrt(M_p*M_p + p_p*p_p);
    G4double KEp = Ep - M_p;

    if (KEp > 10*MeV && p_p > 0 && PassesPauliBlocking(p_p)) {
        fParticleGun->SetParticleDefinition(G4Proton::Proton());
        fParticleGun->SetParticleEnergy(KEp);
        fParticleGun->SetParticleMomentumDirection(p_p_vec.unit());
        fParticleGun->SetParticlePosition(position);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    eventInfo->SetPrimaryVertex(position);
    eventInfo->SetPrimaryMomentum(p_e_vec/GeV);
    eventInfo->SetPrimaryEnergy(Ee / GeV);
    eventInfo->SetPrimaryPDG(11);

    G4double W2_val = M_n*M_n + 2*M_n*(Enu-Ee) - Q2;
    G4double W = (W2_val > 0.01*GeV*GeV) ? sqrt(W2_val) : 0.939*GeV;

    eventInfo->SetNeutrinoKinematics(Q2/(GeV*GeV), W/GeV, Enu/GeV, "NuE_CCQE");
}

void PrimaryGeneratorAction::GenerateNueCC1Pi(G4Event* anEvent, EventInformation* eventInfo, G4double Enu)
{
    if (Enu < 0.30*GeV) {
        GenerateNueCCQE(anEvent, eventInfo, Enu);
        return;
    }

    G4ThreeVector position = SampleVertexInTPC();

    G4double m_e = 0.000511*GeV;
    G4double M_N = 0.939*GeV;
    G4double W = SampleResonance_W();

    G4double Q2_scale = 0.5*GeV*GeV;
    G4double Q2 = -Q2_scale * log(G4UniformRand());
    G4double Q2_max = 2.0*M_N*Enu;
    if (Q2 > Q2_max) Q2 = Q2_max * G4UniformRand();
    if (Q2 < 0.02*GeV*GeV) Q2 = 0.02*GeV*GeV;

    G4double nu = (W*W - M_N*M_N + Q2) / (2.0*M_N);
    G4double Ee = Enu - nu;

    if (Ee < m_e + 20*MeV) Ee = m_e + 20*MeV;
    if (Ee > Enu - 200*MeV) Ee = Enu - 200*MeV;

    G4double cos_theta = 1.0 - (Q2 + m_e*m_e)/(2.0*Enu*Ee);
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

    G4double theta_e = acos(cos_theta);
    G4double phi_e = G4UniformRand() * 2.0 * pi;

    G4ThreeVector e_dir(sin(theta_e)*cos(phi_e),
                        sin(theta_e)*sin(phi_e),
                        cos(theta_e));

    G4double KEe = Ee - m_e;
    if (KEe <= 0) KEe = 20*MeV;

    fParticleGun->SetParticleDefinition(G4Electron::Electron());
    fParticleGun->SetParticleEnergy(KEe);
    fParticleGun->SetParticleMomentumDirection(e_dir);
    fParticleGun->SetParticlePosition(position);
    fParticleGun->GeneratePrimaryVertex(anEvent);

    G4double p_e = (Ee*Ee - m_e*m_e > 0) ? sqrt(Ee*Ee - m_e*m_e) : 0.01*GeV;
    G4ThreeVector p_e_vec = e_dir * p_e;

    G4double E_available = Enu - Ee;
    G4double E_remaining = E_available;

    G4double KEp = E_remaining * (0.4 + 0.2*G4UniformRand());
    E_remaining -= KEp;

    G4double KEpi = E_remaining * (0.6 + 0.3*G4UniformRand());

    G4int pionPDG;
    G4double randPion = G4UniformRand();
    if (randPion < 0.55) pionPDG = 211;
    else                 pionPDG = 111;

    G4ThreeVector p_nu_vec_nue(0, 0, Enu);
    G4ThreeVector q_vec_nue = p_nu_vec_nue - p_e_vec;
    G4ThreeVector q_dir_nue = (q_vec_nue.mag() > 0) ? q_vec_nue.unit()
                                                    : G4ThreeVector(0,0,1);

    std::vector<G4int> hadPDGs;
    std::vector<G4double> hadKEs;
    std::vector<G4ThreeVector> hadDirs;

    if (KEp > 20*MeV) {
        hadPDGs.push_back(2212);
        hadKEs.push_back(KEp);
        hadDirs.push_back(SampleForwardDirection(q_dir_nue, 0.35));
    }
    if (KEpi > 30*MeV) {
        hadPDGs.push_back(pionPDG);
        hadKEs.push_back(KEpi);
        hadDirs.push_back(SampleForwardDirection(q_dir_nue, 0.45));
    }

    ApplyFSI(anEvent, position, hadPDGs, hadKEs, hadDirs);

    eventInfo->SetPrimaryVertex(position);
    eventInfo->SetPrimaryMomentum(p_e_vec/GeV);
    eventInfo->SetPrimaryEnergy(Ee / GeV);
    eventInfo->SetPrimaryPDG(11);
    eventInfo->SetNeutrinoKinematics(Q2/(GeV*GeV), W/GeV, Enu/GeV, "NuE_CC1Pi");
}

void PrimaryGeneratorAction::GenerateNueDIS(G4Event* anEvent, EventInformation* eventInfo, G4double Enu)
{
    if (Enu < 0.5*GeV) {
        GenerateNueCCQE(anEvent, eventInfo, Enu);
        return;
    }

    G4ThreeVector position = SampleVertexInTPC();

    G4double m_e = 0.000511*GeV;
    G4double M_N = 0.939*GeV;

    G4double y = SampleDIS_y(Enu);
    G4double nu = y * Enu;
    G4double Ee = Enu - nu;

    if (Ee < m_e + 50*MeV) {
        Ee = m_e + 50*MeV;
        nu = Enu - Ee;
    }

    G4double Q2;
    if (G4UniformRand() < 0.7)
        Q2 = -2.0*GeV*GeV * log(G4UniformRand());
    else
        Q2 = 0.5*GeV*GeV + 2.5*GeV*GeV * G4UniformRand();

    if (Q2 > 2.0*M_N*Enu) Q2 = 2.0*M_N*Enu * 0.8;
    if (Q2 < 0.1*GeV*GeV) Q2 = 0.1*GeV*GeV;

    G4double W2 = M_N*M_N + 2.0*M_N*nu - Q2;
    G4double W = (W2 > 0) ? sqrt(W2) : 1.1*GeV + 0.3*GeV * G4UniformRand();
    if (W < 1.1*GeV) W = 1.1*GeV + 0.5*GeV * G4UniformRand();

    G4double cos_theta = 1.0 - (Q2 + m_e*m_e)/(2.0*Enu*Ee);
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

    G4double theta_e = acos(cos_theta);
    G4double phi_e = G4UniformRand() * 2.0 * pi;

    G4ThreeVector e_dir(sin(theta_e)*cos(phi_e),
                        sin(theta_e)*sin(phi_e),
                        cos(theta_e));

    G4double KEe = Ee - m_e;
    if (KEe <= 0) KEe = 50*MeV;

    fParticleGun->SetParticleDefinition(G4Electron::Electron());
    fParticleGun->SetParticleEnergy(KEe);
    fParticleGun->SetParticleMomentumDirection(e_dir);
    fParticleGun->SetParticlePosition(position);
    fParticleGun->GeneratePrimaryVertex(anEvent);

    G4double p_e = (Ee*Ee - m_e*m_e > 0) ? sqrt(Ee*Ee - m_e*m_e) : 0.01*GeV;
    G4ThreeVector p_e_vec = e_dir * p_e;

    G4double E_available = Enu - Ee;
    G4int nParticles = 5 + G4int(G4UniformRand() * 5);

    std::vector<G4ParticleDefinition*> hadrons = {
        G4Proton::Proton(), G4Neutron::Neutron(),
        G4PionPlus::PionPlus(), G4PionMinus::PionMinus(),
        G4PionZero::PionZero()
    };

    G4ThreeVector p_nu_vec_dis(0, 0, Enu);
    G4ThreeVector q_vec = p_nu_vec_dis - p_e_vec;
    G4ThreeVector q_dir = (q_vec.mag() > 0) ? q_vec.unit() : G4ThreeVector(0,0,1);

    G4double E_per = E_available / nParticles;
    G4double E_rem = E_available;

    for (G4int i = 0; i < nParticles; i++) {
        G4double frac = exp(-2.5 * G4UniformRand());
        G4double E_had = E_per * frac * 2.5;
        if (E_had < 70*MeV) E_had = 70*MeV;
        if (E_had > E_rem * 0.7) E_had = E_rem * 0.7;
        if (i == nParticles - 1) E_had = E_rem * 0.8;
        E_rem -= E_had;
        if (E_rem < 0) E_rem = 0;

        if (E_had > 60*MeV) {
            G4ParticleDefinition* hadron;
            G4double rH = G4UniformRand();
            if      (rH < 0.22) hadron = G4Proton::Proton();
            else if (rH < 0.44) hadron = G4Neutron::Neutron();
            else if (rH < 0.63) hadron = G4PionPlus::PionPlus();
            else if (rH < 0.82) hadron = G4PionMinus::PionMinus();
            else                hadron = G4PionZero::PionZero();

            G4double thetaScale = (i == 0) ? 0.25 : 0.40;
            G4ThreeVector had_dir = SampleForwardDirection(q_dir, thetaScale);

            fParticleGun->SetParticleDefinition(hadron);
            fParticleGun->SetParticleEnergy(E_had);
            fParticleGun->SetParticleMomentumDirection(had_dir);
            fParticleGun->SetParticlePosition(position);
            fParticleGun->GeneratePrimaryVertex(anEvent);
        }
    }

    eventInfo->SetPrimaryVertex(position);
    eventInfo->SetPrimaryMomentum(p_e_vec/GeV);
    eventInfo->SetPrimaryEnergy(Ee / GeV);
    eventInfo->SetPrimaryPDG(11);
    eventInfo->SetNeutrinoKinematics(Q2/(GeV*GeV), W/GeV, Enu/GeV, "NuE_DIS");
}

void PrimaryGeneratorAction::GenerateNuTauCCQE(G4Event* anEvent, EventInformation* eventInfo, G4double Enu)
{
    G4ThreeVector position = SampleVertexInTPC();

    G4double m_tau = 1.777*GeV;
    G4double M_n = 0.939*GeV;
    G4double M_p = 0.938*GeV;

    G4double Q2 = SampleCCQE_Q2(Enu);

    G4double p_fermi = 0.220*GeV;
    G4ThreeVector pF_vec = G4RandomDirection() * p_fermi;
    G4double Eb = 27*MeV;

    G4double Etau = Enu - Q2/(2.0*M_n);

    if (Etau < m_tau + 50*MeV) Etau = m_tau + 50*MeV;
    if (Etau > Enu - 100*MeV) Etau = Enu - 100*MeV;

    G4double cos_theta = 1.0 - (Q2 + m_tau*m_tau)/(2.0*Enu*Etau);
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

    G4double theta_tau = acos(cos_theta);
    G4double phi_tau = G4UniformRand() * 2.0 * pi;

    G4ThreeVector tau_dir(sin(theta_tau)*cos(phi_tau),
                          sin(theta_tau)*sin(phi_tau),
                          cos(theta_tau));

    SimulateTauDecay(anEvent, position, Etau, tau_dir);

    G4double p_tau = (Etau*Etau - m_tau*m_tau > 0) ? sqrt(Etau*Etau - m_tau*m_tau) : 0.1*GeV;
    G4ThreeVector p_tau_vec = tau_dir * p_tau;
    G4ThreeVector p_nu_vec(0, 0, Enu);
    G4ThreeVector p_p_vec = p_nu_vec + pF_vec - p_tau_vec;

    G4double p_p = p_p_vec.mag();
    G4double KEp = sqrt(M_p*M_p + p_p*p_p) - M_p;

    if (KEp > 10*MeV && p_p > 0) {
        fParticleGun->SetParticleDefinition(G4Proton::Proton());
        fParticleGun->SetParticleEnergy(KEp);
        fParticleGun->SetParticleMomentumDirection(p_p_vec.unit());
        fParticleGun->SetParticlePosition(position);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    eventInfo->SetPrimaryVertex(position);
    eventInfo->SetPrimaryMomentum(p_tau_vec/GeV);
    eventInfo->SetPrimaryEnergy(Etau / GeV);
    eventInfo->SetPrimaryPDG(15);

    G4double W2_val = M_n*M_n + 2*M_n*(Enu-Etau) - Q2;
    G4double W = (W2_val > 0.01*GeV*GeV) ? sqrt(W2_val) : M_n;

    eventInfo->SetNeutrinoKinematics(Q2/(GeV*GeV), W/GeV, Enu/GeV, "NuTau_CCQE");
}

void PrimaryGeneratorAction::GenerateNuTauDIS(G4Event* anEvent, EventInformation* eventInfo, G4double Enu)
{
    G4ThreeVector position = SampleVertexInTPC();

    G4double m_tau = 1.777*GeV;
    G4double M_N = 0.939*GeV;

    G4double y = SampleDIS_y(Enu);
    G4double nu = y * Enu;
    G4double Etau = Enu - nu;

    if (Etau < m_tau + 200*MeV) {
        Etau = m_tau + 200*MeV;
        nu = Enu - Etau;
    }

    G4double Q2 = -2.0*GeV*GeV * log(G4UniformRand());
    if (Q2 > 2.0*M_N*Enu) Q2 = 2.0*M_N*Enu * 0.8;
    if (Q2 < 0.1*GeV*GeV) Q2 = 0.1*GeV*GeV;

    G4double W2 = M_N*M_N + 2.0*M_N*nu - Q2;
    G4double W = (W2 > 0) ? sqrt(W2) : 1.5*GeV;

    G4double cos_theta = 1.0 - (Q2 + m_tau*m_tau)/(2.0*Enu*Etau);
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

    G4double theta_tau = acos(cos_theta);
    G4double phi_tau = G4UniformRand() * 2.0 * pi;

    G4ThreeVector tau_dir(sin(theta_tau)*cos(phi_tau),
                          sin(theta_tau)*sin(phi_tau),
                          cos(theta_tau));

    SimulateTauDecay(anEvent, position, Etau, tau_dir);

    G4double p_tau = (Etau*Etau - m_tau*m_tau > 0) ? sqrt(Etau*Etau - m_tau*m_tau) : 0.1*GeV;
    G4ThreeVector p_tau_vec = tau_dir * p_tau;

    G4double E_available = nu;
    G4int nParticles = 3 + G4int(G4UniformRand() * 4);

    std::vector<G4ParticleDefinition*> hadrons = {
        G4Proton::Proton(), G4Neutron::Neutron(),
        G4PionPlus::PionPlus(), G4PionMinus::PionMinus(),
        G4PionZero::PionZero()
    };

    G4double E_rem = E_available;
    for (G4int i = 0; i < nParticles && E_rem > 100*MeV; i++) {
        G4double E_had = E_rem * (0.2 + 0.3*G4UniformRand());
        if (i == nParticles - 1) E_had = E_rem * 0.8;
        E_rem -= E_had;

        if (E_had > 60*MeV) {
            G4ParticleDefinition* hadron = hadrons[G4int(G4UniformRand() * hadrons.size())];
            fParticleGun->SetParticleDefinition(hadron);
            fParticleGun->SetParticleEnergy(E_had);
            fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
            fParticleGun->SetParticlePosition(position);
            fParticleGun->GeneratePrimaryVertex(anEvent);
        }
    }

    eventInfo->SetPrimaryVertex(position);
    eventInfo->SetPrimaryMomentum(p_tau_vec/GeV);
    eventInfo->SetPrimaryEnergy(Etau / GeV);
    eventInfo->SetPrimaryPDG(15);
    eventInfo->SetNeutrinoKinematics(Q2/(GeV*GeV), W/GeV, Enu/GeV, "NuTau_DIS");
}

void PrimaryGeneratorAction::SimulateTauDecay(G4Event* anEvent,
                                               const G4ThreeVector& position,
                                               G4double tauEnergy,
                                               const G4ThreeVector& tauDirection)
{
    G4double m_tau = 1.777*GeV;
    G4double p_tau = (tauEnergy*tauEnergy - m_tau*m_tau > 0) ?
                     sqrt(tauEnergy*tauEnergy - m_tau*m_tau) : 0.1*GeV;

    G4double gamma = tauEnergy / m_tau;
    G4double beta = p_tau / tauEnergy;

    G4double rand = G4UniformRand();

    if (rand < 0.178) {
        G4double Ee = tauEnergy * (0.15 + 0.35 * G4UniformRand());
        G4double m_e = 0.000511*GeV;
        G4double KEe = Ee - m_e;
        if (KEe < 10*MeV) KEe = 10*MeV;

        G4double spreadAngle = 0.3 / gamma;
        G4double theta_spread = spreadAngle * sqrt(-2.0 * log(G4UniformRand()));
        G4double phi_spread = G4UniformRand() * 2.0 * pi;

        G4ThreeVector e_dir = tauDirection;
        e_dir.rotateY(theta_spread);
        e_dir.rotateZ(phi_spread);

        fParticleGun->SetParticleDefinition(G4Electron::Electron());
        fParticleGun->SetParticleEnergy(KEe);
        fParticleGun->SetParticleMomentumDirection(e_dir);
        fParticleGun->SetParticlePosition(position);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
    else if (rand < 0.352) {
        G4double Emu = tauEnergy * (0.15 + 0.35 * G4UniformRand());
        G4double m_mu = 0.1057*GeV;
        G4double KEmu = Emu - m_mu;
        if (KEmu < 50*MeV) KEmu = 50*MeV;

        G4double spreadAngle = 0.3 / gamma;
        G4double theta_spread = spreadAngle * sqrt(-2.0 * log(G4UniformRand()));
        G4double phi_spread = G4UniformRand() * 2.0 * pi;

        G4ThreeVector mu_dir = tauDirection;
        mu_dir.rotateY(theta_spread);
        mu_dir.rotateZ(phi_spread);

        fParticleGun->SetParticleDefinition(G4MuonMinus::MuonMinus());
        fParticleGun->SetParticleEnergy(KEmu);
        fParticleGun->SetParticleMomentumDirection(mu_dir);
        fParticleGun->SetParticlePosition(position);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
    else if (rand < 0.460) {
        G4double Epi = tauEnergy * (0.3 + 0.4 * G4UniformRand());
        G4double m_pi = 0.1396*GeV;
        G4double KEpi = Epi - m_pi;
        if (KEpi < 50*MeV) KEpi = 50*MeV;

        G4ThreeVector pi_dir = tauDirection;
        G4double theta_spread = 0.2 / gamma * sqrt(-2.0 * log(G4UniformRand()));
        pi_dir.rotateY(theta_spread);
        pi_dir.rotateZ(G4UniformRand() * 2.0 * pi);

        fParticleGun->SetParticleDefinition(G4PionMinus::PionMinus());
        fParticleGun->SetParticleEnergy(KEpi);
        fParticleGun->SetParticleMomentumDirection(pi_dir);
        fParticleGun->SetParticlePosition(position);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
    else if (rand < 0.715) {
        G4double E_total = tauEnergy * (0.4 + 0.3 * G4UniformRand());

        G4double KEpi_minus = E_total * 0.5;
        if (KEpi_minus < 50*MeV) KEpi_minus = 50*MeV;

        G4ThreeVector pi_dir = tauDirection;
        pi_dir.rotateY(0.3/gamma * sqrt(-2.0*log(G4UniformRand())));
        pi_dir.rotateZ(G4UniformRand() * 2.0 * pi);

        fParticleGun->SetParticleDefinition(G4PionMinus::PionMinus());
        fParticleGun->SetParticleEnergy(KEpi_minus);
        fParticleGun->SetParticleMomentumDirection(pi_dir);
        fParticleGun->SetParticlePosition(position);
        fParticleGun->GeneratePrimaryVertex(anEvent);

        G4double KEpi_zero = E_total * 0.4;
        if (KEpi_zero < 50*MeV) KEpi_zero = 50*MeV;

        G4ThreeVector pi0_dir = tauDirection;
        pi0_dir.rotateY(0.4/gamma * sqrt(-2.0*log(G4UniformRand())));
        pi0_dir.rotateZ(G4UniformRand() * 2.0 * pi);

        fParticleGun->SetParticleDefinition(G4PionZero::PionZero());
        fParticleGun->SetParticleEnergy(KEpi_zero);
        fParticleGun->SetParticleMomentumDirection(pi0_dir);
        fParticleGun->SetParticlePosition(position);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
    else {
        G4double E_total = tauEnergy * (0.4 + 0.3 * G4UniformRand());
        G4int nPions = 3;

        G4ParticleDefinition* pions[] = {
            G4PionMinus::PionMinus(),
            G4PionPlus::PionPlus(),
            G4PionMinus::PionMinus()
        };

        G4double E_per = E_total / nPions;
        for (G4int i = 0; i < nPions; i++) {
            G4double KEpi = E_per * (0.5 + G4UniformRand());
            if (KEpi < 40*MeV) KEpi = 40*MeV;

            G4ThreeVector dir = tauDirection;
            dir.rotateY(0.4/gamma * sqrt(-2.0*log(G4UniformRand())));
            dir.rotateZ(G4UniformRand() * 2.0 * pi);

            fParticleGun->SetParticleDefinition(pions[i]);
            fParticleGun->SetParticleEnergy(KEpi);
            fParticleGun->SetParticleMomentumDirection(dir);
            fParticleGun->SetParticlePosition(position);
            fParticleGun->GeneratePrimaryVertex(anEvent);
        }
    }
}

void PrimaryGeneratorAction::GenerateNCQE(G4Event* anEvent, EventInformation* eventInfo, G4double Enu)
{
    G4ThreeVector position = SampleVertexInTPC();

    G4double M_N = 0.939*GeV;

    G4double Q2 = SampleCCQE_Q2(Enu) * 0.7;

    G4double p_fermi = 0.220*GeV;
    G4ThreeVector pF_vec = G4RandomDirection() * p_fermi;

    G4double T_recoil = Q2 / (2.0 * M_N);

    G4double rand = G4UniformRand();
    if (rand < 0.6) {
        if (T_recoil > 20*MeV) {
            G4ThreeVector recoil_dir = G4RandomDirection();
            G4ThreeVector q_dir(0, 0, 1);
            recoil_dir = (q_dir + 0.5 * G4RandomDirection()).unit();

            fParticleGun->SetParticleDefinition(G4Proton::Proton());
            fParticleGun->SetParticleEnergy(T_recoil);
            fParticleGun->SetParticleMomentumDirection(recoil_dir);
            fParticleGun->SetParticlePosition(position);
            fParticleGun->GeneratePrimaryVertex(anEvent);
        }
    }
    if (G4UniformRand() < 0.15 && Enu > 0.5*GeV) {
        G4double KEpi = 50*MeV + 100*MeV * G4UniformRand();
        fParticleGun->SetParticleDefinition(G4PionPlus::PionPlus());
        fParticleGun->SetParticleEnergy(KEpi);
        fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
        fParticleGun->SetParticlePosition(position);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    eventInfo->SetPrimaryVertex(position);
    eventInfo->SetPrimaryMomentum(G4ThreeVector(0, 0, 0));
    eventInfo->SetPrimaryEnergy(0.0);
    eventInfo->SetPrimaryPDG(0);

    G4double W = M_N;
    eventInfo->SetNeutrinoKinematics(Q2/(GeV*GeV), W/GeV, Enu/GeV, "NC_QE");
}

void PrimaryGeneratorAction::GenerateNCRes(G4Event* anEvent, EventInformation* eventInfo, G4double Enu)
{
    if (Enu < 0.35*GeV) {
        GenerateNCQE(anEvent, eventInfo, Enu);
        return;
    }

    G4ThreeVector position = SampleVertexInTPC();
    G4double M_N = 0.939*GeV;

    G4double W = SampleResonance_W();
    G4double Q2 = 0.3*GeV*GeV * (-log(G4UniformRand()));
    if (Q2 > 2.0*M_N*Enu) Q2 = M_N*Enu * G4UniformRand();

    G4double E_hadronic = Enu * (0.2 + 0.3 * G4UniformRand());

    G4double KEp = E_hadronic * (0.4 + 0.2 * G4UniformRand());
    G4double KEpi = E_hadronic - KEp;
    if (KEpi < 30*MeV) KEpi = 30*MeV;

    G4int pionPDG;
    G4double randPion = G4UniformRand();
    if (randPion < 0.50)      pionPDG = 111;
    else if (randPion < 0.75) pionPDG = 211;
    else                      pionPDG = -211;

    std::vector<G4int> hadPDGs;
    std::vector<G4double> hadKEs;
    std::vector<G4ThreeVector> hadDirs;

    if (KEp > 20*MeV) {
        hadPDGs.push_back(2212);
        hadKEs.push_back(KEp);
        hadDirs.push_back(G4RandomDirection());
    }
    if (KEpi > 30*MeV) {
        hadPDGs.push_back(pionPDG);
        hadKEs.push_back(KEpi);
        hadDirs.push_back(G4RandomDirection());
    }

    ApplyFSI(anEvent, position, hadPDGs, hadKEs, hadDirs);

    eventInfo->SetPrimaryVertex(position);
    eventInfo->SetPrimaryMomentum(G4ThreeVector(0, 0, 0));
    eventInfo->SetPrimaryEnergy(0.0);
    eventInfo->SetPrimaryPDG(0);
    eventInfo->SetNeutrinoKinematics(Q2/(GeV*GeV), W/GeV, Enu/GeV, "NC_Res");
}

void PrimaryGeneratorAction::GenerateNCDIS(G4Event* anEvent, EventInformation* eventInfo, G4double Enu)
{
    if (Enu < 0.8*GeV) {
        GenerateNCRes(anEvent, eventInfo, Enu);
        return;
    }

    G4ThreeVector position = SampleVertexInTPC();
    G4double M_N = 0.939*GeV;

    G4double y = SampleDIS_y(Enu);
    y = y * 0.6;

    G4double nu = y * Enu;

    G4double Q2 = -1.5*GeV*GeV * log(G4UniformRand());
    if (Q2 > 2.0*M_N*nu) Q2 = M_N*nu;
    if (Q2 < 0.05*GeV*GeV) Q2 = 0.05*GeV*GeV;

    G4double W2 = M_N*M_N + 2.0*M_N*nu - Q2;
    G4double W = (W2 > 0) ? sqrt(W2) : 1.5*GeV;

    G4double E_hadronic = nu;
    G4int nParticles = 3 + G4int(G4UniformRand() * 4);

    std::vector<G4ParticleDefinition*> hadrons = {
        G4Proton::Proton(), G4Neutron::Neutron(),
        G4PionPlus::PionPlus(), G4PionMinus::PionMinus(),
        G4PionZero::PionZero()
    };

    G4ThreeVector beam_dir(0, 0, 1);

    G4double E_rem = E_hadronic;
    for (G4int i = 0; i < nParticles && E_rem > 80*MeV; i++) {
        G4double E_had = E_rem * (0.2 + 0.3*G4UniformRand());
        if (i == nParticles - 1) E_had = E_rem * 0.8;
        E_rem -= E_had;

        if (E_had > 50*MeV) {
            G4ParticleDefinition* hadron;
            G4double rH = G4UniformRand();
            if      (rH < 0.22) hadron = G4Proton::Proton();
            else if (rH < 0.44) hadron = G4Neutron::Neutron();
            else if (rH < 0.63) hadron = G4PionPlus::PionPlus();
            else if (rH < 0.82) hadron = G4PionMinus::PionMinus();
            else                hadron = G4PionZero::PionZero();

            G4double thetaScale = (i == 0) ? 0.35 : 0.55;
            G4ThreeVector had_dir = SampleForwardDirection(beam_dir, thetaScale);

            fParticleGun->SetParticleDefinition(hadron);
            fParticleGun->SetParticleEnergy(E_had);
            fParticleGun->SetParticleMomentumDirection(had_dir);
            fParticleGun->SetParticlePosition(position);
            fParticleGun->GeneratePrimaryVertex(anEvent);
        }
    }

    eventInfo->SetPrimaryVertex(position);
    eventInfo->SetPrimaryMomentum(G4ThreeVector(0, 0, 0));
    eventInfo->SetPrimaryEnergy(0.0);
    eventInfo->SetPrimaryPDG(0);
    eventInfo->SetNeutrinoKinematics(Q2/(GeV*GeV), W/GeV, Enu/GeV, "NC_DIS");
}

void PrimaryGeneratorAction::Generate2p2hMEC(G4Event* anEvent,
                                              EventInformation* eventInfo,
                                              G4double Enu, G4int leptonPDG)
{
    G4ThreeVector position = SampleVertexInTPC();

    G4double m_l;
    G4ParticleDefinition* lepton;
    if (leptonPDG == 13) {
        m_l = 0.1057*GeV;
        lepton = G4MuonMinus::MuonMinus();
    } else {
        m_l = 0.000511*GeV;
        lepton = G4Electron::Electron();
    }

    G4double M_N = 0.939*GeV;

    G4double omega = 0.0;
    G4int tries = 0;
    while (tries < 100) {
        omega = std::abs(G4RandGauss::shoot(0.12*GeV, 0.08*GeV));
        if (omega > 0.03*GeV && omega < 0.6*GeV && omega < Enu - m_l) break;
        tries++;
    }
    if (omega < 0.03*GeV) omega = 0.05*GeV;

    G4double El = Enu - omega;
    if (El < m_l + 10*MeV) El = m_l + 10*MeV;

    G4double Q2 = 0.3*GeV*GeV * (-std::log(G4UniformRand() + 1e-10));
    G4double Q2_max = 2.0 * Enu * El;
    if (Q2 > Q2_max) Q2 = Q2_max * G4UniformRand();
    if (Q2 < 0.01*GeV*GeV) Q2 = 0.01*GeV*GeV;

    G4double cos_theta = 1.0 - (Q2 + m_l*m_l) / (2.0*Enu*El);
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

    G4double theta_l = std::acos(cos_theta);
    G4double phi_l = G4UniformRand() * 2.0 * pi;

    G4ThreeVector l_dir(std::sin(theta_l)*std::cos(phi_l),
                        std::sin(theta_l)*std::sin(phi_l),
                        std::cos(theta_l));

    G4double KEl = El - m_l;
    if (KEl <= 0) KEl = 10*MeV;

    fParticleGun->SetParticleDefinition(lepton);
    fParticleGun->SetParticleEnergy(KEl);
    fParticleGun->SetParticleMomentumDirection(l_dir);
    fParticleGun->SetParticlePosition(position);
    fParticleGun->GeneratePrimaryVertex(anEvent);

    G4double E_available = omega;
    G4double Eb = 27*MeV;
    E_available -= 2.0 * Eb;
    if (E_available < 20*MeV) E_available = 20*MeV;

    G4double splitFrac = 0.3 + 0.4 * G4UniformRand();
    G4double KE_p1 = E_available * splitFrac;
    G4double KE_p2 = E_available * (1.0 - splitFrac);

    G4ThreeVector q_dir(0, 0, 1);
    G4double openAngle = (60.0 + 60.0 * G4UniformRand()) * M_PI / 180.0;

    G4double theta_p1 = openAngle * 0.5 * (0.5 + 0.5 * G4UniformRand());
    G4double phi_p1 = G4UniformRand() * 2.0 * pi;
    G4ThreeVector p1_dir = q_dir;
    p1_dir.rotateY(theta_p1);
    p1_dir.rotateZ(phi_p1);

    if (KE_p1 > 15*MeV) {
        fParticleGun->SetParticleDefinition(G4Proton::Proton());
        fParticleGun->SetParticleEnergy(KE_p1);
        fParticleGun->SetParticleMomentumDirection(p1_dir);
        fParticleGun->SetParticlePosition(position);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    G4double theta_p2 = openAngle * 0.5 * (0.5 + 0.5 * G4UniformRand());
    G4ThreeVector p2_dir = q_dir;
    p2_dir.rotateY(-theta_p2);
    p2_dir.rotateZ(phi_p1 + M_PI + 0.3*(G4UniformRand()-0.5));

    G4ParticleDefinition* nucleon2 = (G4UniformRand() < 0.80) ?
        static_cast<G4ParticleDefinition*>(G4Proton::Proton()) :
        static_cast<G4ParticleDefinition*>(G4Neutron::Neutron());

    if (KE_p2 > 15*MeV) {
        fParticleGun->SetParticleDefinition(nucleon2);
        fParticleGun->SetParticleEnergy(KE_p2);
        fParticleGun->SetParticleMomentumDirection(p2_dir);
        fParticleGun->SetParticlePosition(position);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    G4double p_l = (El*El - m_l*m_l > 0) ? std::sqrt(El*El - m_l*m_l) : 0.01*GeV;
    G4ThreeVector p_l_vec = l_dir * p_l;

    eventInfo->SetPrimaryVertex(position);
    eventInfo->SetPrimaryMomentum(p_l_vec/GeV);
    eventInfo->SetPrimaryEnergy(El / GeV);
    eventInfo->SetPrimaryPDG(leptonPDG);

    G4String intType = (leptonPDG == 13) ? "NuMu_2p2h" : "NuE_2p2h";
    eventInfo->SetNeutrinoKinematics(Q2/(GeV*GeV), 0.0, Enu/GeV, intType);
}

G4bool PrimaryGeneratorAction::PassesPauliBlocking(G4double nucleonMomentum) const
{
    const G4double kF = 0.220*GeV;
    return (nucleonMomentum > kF);
}

void PrimaryGeneratorAction::ApplyFSI(G4Event* anEvent,
                                       const G4ThreeVector& vertex,
                                       std::vector<G4int>& hadronPDGs,
                                       std::vector<G4double>& hadronKEs,
                                       std::vector<G4ThreeVector>& hadronDirs)
{
    const G4double R_nucleus = 3.9;
    const G4double rho0 = 0.16;

    std::vector<G4int> newPDGs;
    std::vector<G4double> newKEs;
    std::vector<G4ThreeVector> newDirs;

    std::vector<G4int> knockoutPDGs;
    std::vector<G4double> knockoutKEs;
    std::vector<G4ThreeVector> knockoutDirs;

    for (size_t i = 0; i < hadronPDGs.size(); i++) {
        G4int pdg = hadronPDGs[i];
        G4double KE = hadronKEs[i];
        G4ThreeVector dir = hadronDirs[i];
        G4int absPDG = std::abs(pdg);

        if (absPDG == 211 || absPDG == 111) {
            G4double rand = G4UniformRand();

            G4double T_pi = KE;
            G4double pAbsorb, pCEX, pQE;

            if (T_pi < 50*MeV) {
                pAbsorb = 0.35; pCEX = 0.10; pQE = 0.25;
            } else if (T_pi < 200*MeV) {
                pAbsorb = 0.25; pCEX = 0.12; pQE = 0.30;
            } else if (T_pi < 400*MeV) {
                pAbsorb = 0.15; pCEX = 0.08; pQE = 0.25;
            } else {
                pAbsorb = 0.08; pCEX = 0.05; pQE = 0.15;
            }

            if (rand < pAbsorb) {
                fFSIAbsorptionCount++;
                G4double KE_n1 = KE * (0.3 + 0.4*G4UniformRand());
                G4double KE_n2 = KE - KE_n1;

                if (KE_n1 > 20*MeV) {
                    knockoutPDGs.push_back(2212);
                    knockoutKEs.push_back(KE_n1);
                    knockoutDirs.push_back(G4RandomDirection());
                }
                if (KE_n2 > 20*MeV) {
                    knockoutPDGs.push_back((G4UniformRand() < 0.5) ? 2212 : 2112);
                    knockoutKEs.push_back(KE_n2);
                    knockoutDirs.push_back(G4RandomDirection());
                }
                continue;
            }
            else if (rand < pAbsorb + pCEX) {
                fFSIChargeExchangeCount++;
                if (pdg == 211) {
                    pdg = 111;
                    knockoutPDGs.push_back(2212);
                    knockoutKEs.push_back(30*MeV + 50*MeV*G4UniformRand());
                    knockoutDirs.push_back(G4RandomDirection());
                } else if (pdg == -211) {
                    pdg = 111;
                } else if (pdg == 111) {
                    if (G4UniformRand() < 0.55) {
                        pdg = -211;
                        knockoutPDGs.push_back(2212);
                        knockoutKEs.push_back(20*MeV + 40*MeV*G4UniformRand());
                        knockoutDirs.push_back(G4RandomDirection());
                    } else {
                        pdg = 211;
                    }
                }
                KE *= (0.7 + 0.2*G4UniformRand());
            }
            else if (rand < pAbsorb + pCEX + pQE) {
                KE *= (0.5 + 0.4*G4UniformRand());
                G4double deflection = 0.3 + 0.5*G4UniformRand();
                dir.rotateY(deflection);
                dir.rotateZ(G4UniformRand() * 2.0 * M_PI);
            }

            if (KE > 30*MeV) {
                newPDGs.push_back(pdg);
                newKEs.push_back(KE);
                newDirs.push_back(dir);
            }
        }
        else if (absPDG == 2212) {
            G4double rand = G4UniformRand();

            G4double pRescat = (KE < 100*MeV) ? 0.25 : 0.15;

            if (rand < pRescat) {
                KE *= (0.6 + 0.3*G4UniformRand());
                G4double deflection = 0.2 + 0.6*G4UniformRand();
                dir.rotateY(deflection);
                dir.rotateZ(G4UniformRand() * 2.0 * M_PI);

                if (G4UniformRand() < 0.3 && KE > 50*MeV) {
                    G4double KE_knock = 20*MeV + 40*MeV*G4UniformRand();
                    KE -= KE_knock;
                    knockoutPDGs.push_back((G4UniformRand() < 0.5) ? 2212 : 2112);
                    knockoutKEs.push_back(KE_knock);
                    knockoutDirs.push_back(G4RandomDirection());
                }
            }

            G4double M_p = 0.938*GeV;
            G4double p_mag = std::sqrt(KE*KE + 2.0*KE*M_p);
            if (!PassesPauliBlocking(p_mag)) {
                continue;
            }

            if (KE > 15*MeV) {
                newPDGs.push_back(pdg);
                newKEs.push_back(KE);
                newDirs.push_back(dir);
            }
        }
        else if (absPDG == 2112) {
            G4double M_n = 0.939*GeV;
            G4double p_mag = std::sqrt(KE*KE + 2.0*KE*M_n);
            if (PassesPauliBlocking(p_mag) && KE > 10*MeV) {
                newPDGs.push_back(pdg);
                newKEs.push_back(KE);
                newDirs.push_back(dir);
            }
        }
        else {
            newPDGs.push_back(pdg);
            newKEs.push_back(KE);
            newDirs.push_back(dir);
        }
    }

    for (size_t k = 0; k < knockoutPDGs.size(); k++) {
        newPDGs.push_back(knockoutPDGs[k]);
        newKEs.push_back(knockoutKEs[k]);
        newDirs.push_back(knockoutDirs[k]);
    }

    hadronPDGs = std::move(newPDGs);
    hadronKEs = std::move(newKEs);
    hadronDirs = std::move(newDirs);

    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    for (size_t i = 0; i < hadronPDGs.size(); i++) {
        G4ParticleDefinition* particle = pTable->FindParticle(hadronPDGs[i]);
        if (!particle) continue;
        if (hadronKEs[i] < 10*MeV) continue;

        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticleEnergy(hadronKEs[i]);
        fParticleGun->SetParticleMomentumDirection(hadronDirs[i]);
        fParticleGun->SetParticlePosition(vertex);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
}

void PrimaryGeneratorAction::GenerateCosmicRay(G4Event* anEvent, EventInformation* eventInfo)
{
    G4double energy = SampleCosmicEnergy();
    G4double theta = SampleCosmicAngle();
    G4double phi = G4UniformRand() * 2.0 * pi;

    G4ThreeVector direction(sin(theta)*cos(phi), -cos(theta), sin(theta)*sin(phi));
    G4ThreeVector initialPosition = SampleCosmicPosition();
    G4ThreeVector entryPoint = FindLArEntryPoint(initialPosition, direction);

    G4ParticleDefinition* particle = nullptr;
    G4double randParticle = G4UniformRand();

    if (randParticle < 0.79) {
        particle = (G4UniformRand() < 0.5)
            ? static_cast<G4ParticleDefinition*>(G4MuonMinus::MuonMinus())
            : static_cast<G4ParticleDefinition*>(G4MuonPlus::MuonPlus());
    } else if (randParticle < 0.92) {
        particle = G4Proton::Proton();
    } else if (randParticle < 0.99) {
        particle = G4Neutron::Neutron();
    } else {
        particle = (G4UniformRand() < 0.5)
            ? static_cast<G4ParticleDefinition*>(G4Electron::Electron())
            : static_cast<G4ParticleDefinition*>(G4Positron::Positron());
    }

    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleEnergy(energy);
    fParticleGun->SetParticleMomentumDirection(direction);
    fParticleGun->SetParticlePosition(entryPoint);
    fParticleGun->SetParticleTime(G4UniformRand() * GetMaxVisiblePrimaryTime());
    fParticleGun->GeneratePrimaryVertex(anEvent);

    G4double p_mag = sqrt(energy*energy + 2*energy*particle->GetPDGMass());
    G4ThreeVector p_vec = direction * p_mag;

    eventInfo->SetPrimaryVertex(entryPoint);
    eventInfo->SetPrimaryMomentum(p_vec/GeV);
    eventInfo->SetPrimaryEnergy((energy + particle->GetPDGMass())/GeV);
    eventInfo->SetPrimaryPDG(particle->GetPDGEncoding());
    eventInfo->SetNeutrinoKinematics(-1, -1, -1, "Cosmic");
}

void PrimaryGeneratorAction::GenerateCosmicOverlay(G4Event* anEvent)
{
    const G4double maxPrimaryTime = GetMaxVisiblePrimaryTime();
    G4int nOverlay = static_cast<G4int>(std::round(CLHEP::RandPoisson::shoot(fCosmicOverlayMean)));
    nOverlay = std::max(0, std::min(nOverlay, 60));

    for (G4int i = 0; i < nOverlay; i++) {
        G4double energy = SampleCosmicEnergy();
        G4double theta = SampleCosmicAngle();
        G4double phi = G4UniformRand() * 2.0 * pi;

        G4ThreeVector direction(sin(theta)*cos(phi), -cos(theta), sin(theta)*sin(phi));
        G4ThreeVector initialPosition = SampleCosmicPosition();
        G4ThreeVector entryPoint = FindLArEntryPoint(initialPosition, direction);

        G4ParticleDefinition* particle = nullptr;
        G4double randParticle = G4UniformRand();
        if (randParticle < 0.79) {
            particle = (G4UniformRand() < 0.5)
                ? static_cast<G4ParticleDefinition*>(G4MuonMinus::MuonMinus())
                : static_cast<G4ParticleDefinition*>(G4MuonPlus::MuonPlus());
        } else if (randParticle < 0.92) {
            particle = G4Proton::Proton();
        } else if (randParticle < 0.99) {
            particle = G4Neutron::Neutron();
        } else {
            particle = (G4UniformRand() < 0.5)
                ? static_cast<G4ParticleDefinition*>(G4Electron::Electron())
                : static_cast<G4ParticleDefinition*>(G4Positron::Positron());
        }

        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticleEnergy(energy);
        fParticleGun->SetParticleMomentumDirection(direction);
        fParticleGun->SetParticlePosition(entryPoint);
        fParticleGun->SetParticleTime((maxPrimaryTime > 0.0) ? G4UniformRand() * maxPrimaryTime : 0.0);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    fParticleGun->SetParticleTime(kMicroBooNEBeamSpillTime);
}

void PrimaryGeneratorAction::GenerateTestParticle(G4Event* anEvent)
{
    fParticleGun->SetParticleDefinition(G4MuonMinus::MuonMinus());
    fParticleGun->SetParticleEnergy(2.*GeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    fParticleGun->SetParticlePosition(G4ThreeVector(0,0,-2.*m));
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

G4double PrimaryGeneratorAction::SampleCosmicEnergy()
{
    G4double alpha = 2.7;
    G4double rand = G4UniformRand();
    G4double energy = fCosmicMinEnergy * pow(
        (pow(fCosmicMaxEnergy/fCosmicMinEnergy, 1-alpha) - 1) * rand + 1,
        1./(1-alpha)
    );
    if (energy < 0.5*GeV) energy = 0.5*GeV + G4UniformRand() * 0.5*GeV;
    return energy;
}

G4double PrimaryGeneratorAction::SampleCosmicAngle()
{
    G4double cosMinTheta = cos(fCosmicMaxTheta);
    G4double cosMaxTheta = cos(fCosmicMinTheta);
    G4double randCos = cosMinTheta + G4UniformRand() * (cosMaxTheta - cosMinTheta);
    return acos(sqrt(randCos));
}

G4ThreeVector PrimaryGeneratorAction::SampleCosmicPosition()
{
    G4double height = 3.*m;
    G4double spread = 4.*m;
    return G4ThreeVector(
        (G4UniformRand() - 0.5) * spread,
        height,
        (G4UniformRand() - 0.5) * spread
    );
}

G4ThreeVector PrimaryGeneratorAction::FindLArEntryPoint(
    G4ThreeVector position, G4ThreeVector direction)
{
    G4double halfX = DetectorConstruction::fTPC_X / 2.0 * cm;
    G4double halfY = DetectorConstruction::fTPC_Y / 2.0 * cm;
    G4double halfZ = DetectorConstruction::fTPC_Z / 2.0 * cm;

    G4ThreeVector entryPoint = position;
    G4double tMin = 1e10;
    G4bool foundIntersection = false;

    std::vector<std::pair<G4int, G4double>> faces = {
        {0, -halfX}, {0, halfX},
        {1, -halfY}, {1, halfY},
        {2, -halfZ}, {2, halfZ}
    };

    for (const auto& face : faces) {
        G4int axis = face.first;
        G4double boundary = face.second;
        if (std::abs(direction[axis]) < 1e-6) continue;

        G4double t = (boundary - position[axis]) / direction[axis];
        if (t < 0 || t > tMin) continue;

        G4ThreeVector intersection = position + t * direction;
        bool withinBounds = true;
        if (axis != 0 && std::abs(intersection.x()) > halfX) withinBounds = false;
        if (axis != 1 && std::abs(intersection.y()) > halfY) withinBounds = false;
        if (axis != 2 && std::abs(intersection.z()) > halfZ) withinBounds = false;

        if (withinBounds) {
            tMin = t;
            entryPoint = intersection;
            foundIntersection = true;
        }
    }

    if (!foundIntersection) {
        entryPoint = G4ThreeVector(
            (G4UniformRand() - 0.5) * 2.0 * halfX,
            (G4UniformRand() - 0.5) * 2.0 * halfY,
            halfZ
        );
    }

    if (std::abs(entryPoint.x()) > halfX + 0.1*cm ||
        std::abs(entryPoint.y()) > halfY + 0.1*cm ||
        std::abs(entryPoint.z()) > halfZ + 0.1*cm) {
        entryPoint.setZ(halfZ);
        if (std::abs(entryPoint.x()) > halfX) entryPoint.setX(halfX * (entryPoint.x() > 0 ? 1 : -1));
        if (std::abs(entryPoint.y()) > halfY) entryPoint.setY(halfY * (entryPoint.y() > 0 ? 1 : -1));
    }

    return entryPoint;
}

void PrimaryGeneratorAction::SetCosmicEnergyRange(G4double minE, G4double maxE)
{ fCosmicMinEnergy = minE; fCosmicMaxEnergy = maxE; }

void PrimaryGeneratorAction::SetCosmicAngleRange(G4double minTheta, G4double maxTheta)
{ fCosmicMinTheta = minTheta; fCosmicMaxTheta = maxTheta; }

void PrimaryGeneratorAction::SetNeutrinoEnergyRange(G4double minE, G4double maxE)
{ fNeutrinoMinEnergy = minE; fNeutrinoMaxEnergy = maxE; }

void PrimaryGeneratorAction::InitializeBNBFlux()
{
    for (G4int i = 0; i <= fNBNBBins; i++) {
        fBNBEnergyBins[i] = i * 0.05;
    }

    G4double fluxData[60] = {
        3.09,11.9,15.3,18.3,22.7,25.0,26.7,28.0,29.6,30.9,31.6,31.6,31.2,30.9,30.6,
        29.9,28.7,27.5,26.3,24.9,23.6,22.2,20.6,19.2,17.8,16.2,14.7,13.2,11.7,10.2,
        8.85,7.65,6.50,5.48,4.63,3.83,3.18,2.57,2.10,1.70,1.35,1.11,0.911,0.723,
        0.621,0.535,0.463,0.404,0.367,0.333,0.308,0.289,0.278,0.269,0.258,0.240,
        0.237,0.229,0.225,0.209
    };

    G4double totalFlux = 0.0;
    for (G4int i = 0; i < fNBNBBins; i++) {
        fBNBFluxWeights[i] = fluxData[i];
        totalFlux += fluxData[i];
    }

    fBNBCumulativeWeights[0] = fBNBFluxWeights[0] / totalFlux;
    for (G4int i = 1; i < fNBNBBins; i++) {
        fBNBCumulativeWeights[i] = fBNBCumulativeWeights[i-1] + fBNBFluxWeights[i] / totalFlux;
    }
}

G4double PrimaryGeneratorAction::SampleBNBEnergy()
{
    for (G4int attempt = 0; attempt < 5000; attempt++) {
        G4double rand = G4UniformRand();
        G4int binIndex = 0;
        for (G4int i = 0; i < fNBNBBins; i++) {
            if (rand < fBNBCumulativeWeights[i]) { binIndex = i; break; }
        }
        G4double Emin = fBNBEnergyBins[binIndex];
        G4double Emax = fBNBEnergyBins[binIndex + 1];
        G4double energy = (Emin + G4UniformRand() * (Emax - Emin)) * GeV;
        if (energy >= fNeutrinoMinEnergy && energy <= fNeutrinoMaxEnergy) {
            return energy;
        }
    }

    return std::max(fNeutrinoMinEnergy, std::min(0.8 * GeV, fNeutrinoMaxEnergy));
}

G4double PrimaryGeneratorAction::SampleCCQE_Q2(G4double Enu)
{
    const G4double MA = 1.03*GeV;
    const G4double MA2 = MA*MA;
    const G4double MN = 0.939*GeV;

    G4double Q2_max = 2.0 * MN * Enu;
    if (Q2_max > 3.0*GeV*GeV) Q2_max = 3.0*GeV*GeV;
    if (Q2_max < 0.05*GeV*GeV) return 0.05*GeV*GeV;

    G4double Q2 = 0.0;
    G4bool accepted = false;
    G4int tries = 0;

    while (!accepted && tries < 2000) {
        tries++;

        G4double u = G4UniformRand();
        Q2 = MA2 * u / (1.0 - u + 1e-10);

        if (Q2 > Q2_max || Q2 < 0.01*GeV*GeV) continue;

        G4double tau = Q2 / (4.0 * MN * MN);
        G4double dipole = 1.0 / (1.0 + Q2/MA2);
        G4double FA2 = 1.2670 * 1.2670 * dipole * dipole * dipole * dipole;

        G4double MV2 = 0.71*GeV*GeV;
        G4double GEV = 1.0 / ((1.0 + Q2/MV2) * (1.0 + Q2/MV2));
        G4double GMV = 4.706 * GEV;

        G4double A = (tau + 1.0) * FA2
                   - (tau - 1.0) * (GEV*GEV + tau*GMV*GMV)/(1.0+tau);
        G4double B = Q2/(MN*MN) * FA2;

        G4double weight = std::abs(A + B);
        G4double maxWeight = 3.0;

        if (G4UniformRand() * maxWeight < weight) {
            accepted = true;
        }
    }

    if (!accepted) {
        G4double u = G4UniformRand();
        Q2 = MA2 * u / (1.0 - 0.9*u);
        if (Q2 > Q2_max) Q2 = Q2_max * G4UniformRand();
    }

    Q2 = std::max(0.02*GeV*GeV, std::min(Q2, Q2_max));
    return Q2;
}

G4double PrimaryGeneratorAction::SampleResonance_W()
{
    const G4double M_Delta = 1.232*GeV;
    const G4double Gamma_Delta = 0.120*GeV;
    G4double W;
    G4bool accept = false;
    G4int tries = 0;
    while (!accept && tries < 1000) {
        W = M_Delta + (G4UniformRand() - 0.5) * 0.8*GeV;
        if (W < 1.08*GeV || W > 2.0*GeV) { tries++; continue; }
        G4double W2 = W*W;
        G4double M2 = M_Delta*M_Delta;
        G4double MG = M_Delta * Gamma_Delta;
        G4double BW = Gamma_Delta / (pow(W2 - M2, 2) + MG*MG);
        if (G4UniformRand() < BW * MG * MG) accept = true;
        tries++;
    }
    if (tries >= 1000) W = 1.232*GeV;
    return W;
}

G4double PrimaryGeneratorAction::SampleDIS_y(G4double Enu)
{
    G4double y;
    G4bool accept = false;
    G4int tries = 0;
    while (!accept && tries < 1000) {
        y = G4UniformRand();
        G4double weight = (1.0 + pow(1.0 - y, 2)) / 2.0;
        if (G4UniformRand() < weight) accept = true;
        tries++;
    }
    if (tries >= 1000) y = 0.3;
    y = std::max(0.1, std::min(0.9, y));
    return y;
}

G4ThreeVector PrimaryGeneratorAction::SampleForwardDirection(
    const G4ThreeVector& axis, G4double thetaScale)
{
    G4ThreeVector n = axis;
    if (n.mag() < 1e-9) n = G4ThreeVector(0, 0, 1);
    else                n = n.unit();

    G4double u = G4UniformRand();
    G4double theta = -thetaScale * std::log(u + 1e-12);
    if (theta > M_PI * 0.5) theta = M_PI * 0.5;
    G4double phi = G4UniformRand() * 2.0 * M_PI;

    G4double sth = std::sin(theta);
    G4ThreeVector local(sth * std::cos(phi),
                        sth * std::sin(phi),
                        std::cos(theta));

    G4ThreeVector ref = (std::fabs(n.z()) < 0.9)
                          ? G4ThreeVector(0, 0, 1)
                          : G4ThreeVector(1, 0, 0);
    G4ThreeVector u1 = (ref.cross(n)).unit();
    G4ThreeVector u2 = n.cross(u1);

    return (local.x() * u1 + local.y() * u2 + local.z() * n).unit();
}

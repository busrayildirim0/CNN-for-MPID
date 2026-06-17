// Global Application State
const state = {
  selectedParticle: 'electron',
  isSimulating: false,
  detectorGeometry: {
    widthX: 2.56,   // Scaled from 256.35 cm
    heightY: 2.33,  // Scaled from 233.0 cm
    lengthZ: 10.37, // Scaled from 1036.8 cm
    anodeX: -1.28,  // Anode wire plane location
    cathodeX: 1.28  // Cathode location
  },
  three: {
    scene: null,
    camera: null,
    renderer: null,
    controls: null,
    cryostat: null,
    larVolume: null,
    anodeGrid: null,
    trackMeshGroup: null,
    driftParticles: []
  },
  canvases: {
    u: null,
    v: null,
    y: null
  },
  simulationData: {
    hits: [],          
    particleType: '',
    vertex: null
  }
};

document.addEventListener('DOMContentLoaded', () => {
  initTabs();
  initParticleSelector();
  initThreeJS();
  initCanvases();
  
  document.getElementById('btn-simulate').addEventListener('click', runSimulationEvent);
  
  clearReadouts();
  logMessage("Detector ready. Choose a particle type and click 'Simulate Event'.");
});

  // UI NAVIGATION & TABS
function initTabs() {
  const tabButtons = document.querySelectorAll('.chart-tab-btn');
  const tabContents = document.querySelectorAll('.chart-tab-content');

  tabButtons.forEach(btn => {
    btn.addEventListener('click', () => {
      const targetId = btn.getAttribute('data-target');
      
      tabButtons.forEach(b => b.classList.remove('active'));
      tabContents.forEach(c => c.classList.remove('active'));
      
      btn.classList.add('active');
      const targetContent = document.getElementById(targetId);
      if (targetContent) {
        targetContent.classList.add('active');
      }
    });
  });
}

function initParticleSelector() {
  const options = document.querySelectorAll('.particle-option');
  options.forEach(option => {
    option.addEventListener('click', () => {
      if (state.isSimulating) return; 
      
      options.forEach(opt => opt.classList.remove('selected'));
      option.classList.add('selected');
      state.selectedParticle = option.getAttribute('data-particle');
      
      logMessage(`Selected particle: ${state.selectedParticle.toUpperCase()}. Ready for simulation.`);
    });
  });
}

function logMessage(text, append = false) {
  const logDiv = document.getElementById('simulation-log');
  if (append) {
    logDiv.innerHTML += `<br>➔ ${text}`;
  } else {
    logDiv.innerHTML = `➔ ${text}`;
  }
  logDiv.scrollTop = logDiv.scrollHeight;
}

   // 2D WIRE CANVAS READOUTS
function initCanvases() {
  const views = ['u', 'v', 'y'];
  views.forEach(v => {
    const canvas = document.getElementById(`canvas-view-${v}`);
    state.canvases[v] = canvas;
    
    canvas.width = 512;
    canvas.height = 512;
    
    clearCanvas(canvas, v.toUpperCase());
  });
}

function clearCanvas(canvas, label) {
  const ctx = canvas.getContext('2d');
  ctx.fillStyle = '#f8fafc'; 
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  
  ctx.strokeStyle = 'rgba(15, 23, 42, 0.04)'; 
  ctx.lineWidth = 1;
  
  const step = 32;
  for (let y = 0; y < canvas.height; y += step) {
    ctx.beginPath();
    ctx.moveTo(0, y);
    ctx.lineTo(canvas.width, y);
    ctx.stroke();
  }
  
  if (label === 'Y') {
    for (let x = 0; x < canvas.width; x += step) {
      ctx.beginPath();
      ctx.moveTo(x, 0);
      ctx.lineTo(x, canvas.height);
      ctx.stroke();
    }
  } else if (label === 'U') {
    // +60 degree wires
    for (let i = -canvas.width; i < canvas.width * 2; i += step) {
      ctx.beginPath();
      ctx.moveTo(i, 0);
      ctx.lineTo(i + canvas.height * 0.577, canvas.height); // tan(30) ~ 0.577
      ctx.stroke();
    }
  } else if (label === 'V') {
    // -60 degree wires
    for (let i = 0; i < canvas.width * 2; i += step) {
      ctx.beginPath();
      ctx.moveTo(i, 0);
      ctx.lineTo(i - canvas.height * 0.577, canvas.height);
      ctx.stroke();
    }
  }
  
  ctx.fillStyle = 'rgba(15, 23, 42, 0.45)'; 
  ctx.font = 'bold 12px Outfit, sans-serif';
  ctx.fillText(`VIEW ${label}`, 15, 25);
  ctx.font = '10px JetBrains Mono, monospace';
  ctx.fillText('Z (Wire #) ➔', 15, canvas.height - 15);
  
  ctx.save();
  ctx.translate(canvas.width - 15, 15);
  ctx.rotate(Math.PI / 2);
  ctx.fillText('X (Time/Drift) ➔', 0, 0);
  ctx.restore();
}

function clearReadouts() {
  const views = ['u', 'v', 'y'];
  views.forEach(v => {
    clearCanvas(state.canvases[v], v.toUpperCase());
  });
}

/**
 * Project a 3D coordinate point (x, y, z) onto U, V, and Y wire views.
 * Anode is at X = -1.28, Cathode at X = 1.28.
 * Ticks correspond to X coordinate (vertical axis in canvas).
 * Wires correspond to Z and Y (horizontal axis in canvas).
 */
function drawProjectedHit(v3, charge, color) {
  const x = v3.x;
  const y = v3.y;
  const z = v3.z;
  
  const widthX = state.detectorGeometry.widthX;
  const heightY = state.detectorGeometry.heightY;
  const lengthZ = state.detectorGeometry.lengthZ;
  
  // 1. Vertical Axis: Time Ticks (maps from X: [-widthX/2, widthX/2] to Canvas Y: [40, 472])
  const tickPerc = (x + widthX/2) / widthX; // [0, 1]
  const canvasY = 40 + tickPerc * 432;
  
  // 2. Horizontal Axis: Wires
  // Y View (Collection wires at 0 degrees, reading along Z coordinate)
  const zPerc = (z + lengthZ/2) / lengthZ; // [0, 1]
  const canvasX_Y = 40 + zPerc * 432;
  
  // U View (+60 degrees wires): Coordinate perpendicular to wires is u = z * cos(30) + y * sin(30)
  const uVal = z * Math.cos(Math.PI / 6) + y * Math.sin(Math.PI / 6);
  const uMax = (lengthZ/2) * Math.cos(Math.PI / 6) + (heightY/2) * Math.sin(Math.PI / 6);
  const uPerc = (uVal + uMax) / (uMax * 2);
  const canvasX_U = 40 + uPerc * 432;
  
  // V View (-60 degrees wires): v = z * cos(30) - y * sin(30)
  const vVal = z * Math.cos(Math.PI / 6) - y * Math.sin(Math.PI / 6);
  const vMax = (lengthZ/2) * Math.cos(Math.PI / 6) + (heightY/2) * Math.sin(Math.PI / 6);
  const vPerc = (vVal + vMax) / (vMax * 2);
  const canvasX_V = 40 + vPerc * 432;
  
  //  glowing dots on respective canvases
  drawGlowingDot(state.canvases.y, canvasX_Y, canvasY, charge, color);
  drawGlowingDot(state.canvases.u, canvasX_U, canvasY, charge, color);
  drawGlowingDot(state.canvases.v, canvasX_V, canvasY, charge, color);
}

function drawGlowingDot(canvas, cx, cy, charge, colorStr) {
  const ctx = canvas.getContext('2d');
  
  const radius = Math.max(2, Math.min(6, charge * 12));
  
  ctx.save();
  const grad = ctx.createRadialGradient(cx, cy, 0, cx, cy, radius * 2);
  grad.addColorStop(0, colorStr);
  grad.addColorStop(0.3, colorStr);
  grad.addColorStop(1, 'rgba(248, 250, 252, 0)');
  
  ctx.fillStyle = grad;
  ctx.beginPath();
  ctx.arc(cx, cy, radius * 2, 0, Math.PI * 2);
  ctx.fill();
  ctx.restore();
}

function initThreeJS() {
  const container = document.getElementById('three-canvas-container');
  const width = container.clientWidth;
  const height = container.clientHeight;
  
  // Scene
  state.three.scene = new THREE.Scene();
  state.three.scene.background = new THREE.Color('#f1f5f9'); // Clean light background
  
  // Camera
  state.three.camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 100);
  state.three.camera.position.set(7, 5, 8); // Perspective showing anode, side, length
  
  // Renderer
  state.three.renderer = new THREE.WebGLRenderer({ antialias: true });
  state.three.renderer.setSize(width, height);
  state.three.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  container.appendChild(state.three.renderer.domElement);
  
  // Controls
  state.three.controls = new THREE.OrbitControls(state.three.camera, state.three.renderer.domElement);
  state.three.controls.enableDamping = true;
  state.three.controls.dampingFactor = 0.05;
  state.three.controls.maxPolarAngle = Math.PI / 2 + 0.1;
  state.three.controls.minDistance = 3;
  state.three.controls.maxDistance = 30;
  
  // Lights
  const ambientLight = new THREE.AmbientLight(0xffffff, 0.45);
  state.three.scene.add(ambientLight);
  
  const dirLight1 = new THREE.DirectionalLight(0xffffff, 0.8);
  dirLight1.position.set(10, 15, 10);
  state.three.scene.add(dirLight1);
  
  const dirLight2 = new THREE.DirectionalLight(0x3b82f6, 0.25);
  dirLight2.position.set(-10, -5, -10);
  state.three.scene.add(dirLight2);

  state.three.trackMeshGroup = new THREE.Group();
  state.three.scene.add(state.three.trackMeshGroup);

  buildDetectorGeometry();
  
  document.getElementById('btn-zoom-in').addEventListener('click', () => {
    state.three.camera.position.multiplyScalar(0.85);
  });
  document.getElementById('btn-zoom-out').addEventListener('click', () => {
    state.three.camera.position.multiplyScalar(1.15);
  });
  document.getElementById('btn-reset-camera').addEventListener('click', () => {
    state.three.camera.position.set(7, 5, 8);
    state.three.controls.target.set(0, 0, 0);
  });
  
  let gridVisible = true;
  const gridHelper = new THREE.GridHelper(20, 20, 0xcb8d9f1, 0xe2e8f0);
  gridHelper.position.y = -state.detectorGeometry.heightY / 2 - 0.05;
  state.three.scene.add(gridHelper);
  
  document.getElementById('btn-toggle-grid').addEventListener('click', () => {
    gridVisible = !gridVisible;
    gridHelper.visible = gridVisible;
  });
  
  let cryoVisible = true;
  document.getElementById('btn-toggle-cryo').addEventListener('click', () => {
    cryoVisible = !cryoVisible;
    state.three.cryostat.visible = cryoVisible;
  });

  window.addEventListener('resize', onWindowResize);
  
  animate();
}

function buildDetectorGeometry() {
  const g = state.detectorGeometry;
  
  const cryoGeo = new THREE.CylinderGeometry(g.heightY * 0.9, g.heightY * 0.9, g.lengthZ + 0.8, 24, 1, true);
  cryoGeo.rotateX(Math.PI / 2);
  
  const cryoMat = new THREE.MeshStandardMaterial({
    color: 0x94a3b8,
    metalness: 0.7,
    roughness: 0.3,
    wireframe: false,
    transparent: true,
    opacity: 0.08,
    side: THREE.DoubleSide
  });
  
  state.three.cryostat = new THREE.Mesh(cryoGeo, cryoMat);
  state.three.scene.add(state.three.cryostat);
  
  const ringGeo = new THREE.TorusGeometry(g.heightY * 0.9, 0.03, 8, 24);
  const ringMat = new THREE.MeshStandardMaterial({ color: 0xcb8d9f1, metalness: 0.8, roughness: 0.2, transparent: true, opacity: 0.4 });
  for (let zOffset = -g.lengthZ/2; zOffset <= g.lengthZ/2; zOffset += g.lengthZ/4) {
    const ring = new THREE.Mesh(ringGeo, ringMat);
    ring.position.z = zOffset;
    state.three.cryostat.add(ring);
  }

  // Liquid Argon Volume (Semi-transparent glowing cyan box)
  const larGeo = new THREE.BoxGeometry(g.widthX, g.heightY, g.lengthZ);
  const larMat = new THREE.MeshStandardMaterial({
    color: 0x0099ff,
    transparent: true,
    opacity: 0.12, 
    roughness: 0.2,
    metalness: 0.1,
    side: THREE.DoubleSide
  });
  state.three.larVolume = new THREE.Mesh(larGeo, larMat);
  state.three.scene.add(state.three.larVolume);
  
  // Glow edge helper for Liquid Argon Volume
  const edges = new THREE.EdgesGeometry(larGeo);
  const lineMat = new THREE.LineBasicMaterial({ color: 0x0891b2, linewidth: 2 });
  const wireframe = new THREE.LineSegments(edges, lineMat);
  state.three.larVolume.add(wireframe);

  // Anode Wire Plane (Grid at X = -1.28)
  const gridHelperYZ = new THREE.GridHelper(g.lengthZ, 50, 0x2563eb, 0xe2e8f0);
  gridHelperYZ.rotateZ(Math.PI / 2); 
  gridHelperYZ.position.x = g.anodeX;
  gridHelperYZ.scale.set(g.heightY / g.lengthZ, 1, 1);
  state.three.scene.add(gridHelperYZ);
  state.three.anodeGrid = gridHelperYZ;
  
  // Cathode Plane (Visual representation at X = +1.28)
  const cathodeGeo = new THREE.PlaneGeometry(g.lengthZ, g.heightY);
  cathodeGeo.rotateY(-Math.PI / 2); // Align parallel to Y-Z
  const cathodeMat = new THREE.MeshBasicMaterial({
    color: 0x3b82f6,
    wireframe: true,
    transparent: true,
    opacity: 0.08
  });
  const cathode = new THREE.Mesh(cathodeGeo, cathodeMat);
  cathode.position.x = g.cathodeX;
  state.three.scene.add(cathode);
  
  // Arrow helper showing Electric Drift Field
  const arrowDir = new THREE.Vector3(-1, 0, 0);
  const arrowOrigin = new THREE.Vector3(1.28, 0, 0);
  const arrowLength = 2.5;
  const arrowColor = 0xd97706; 
  const arrowHelper = new THREE.ArrowHelper(arrowDir, arrowOrigin, arrowLength, arrowColor, 0.3, 0.15);
  state.three.scene.add(arrowHelper);
}

function onWindowResize() {
  const container = document.getElementById('three-canvas-container');
  const width = container.clientWidth;
  const height = container.clientHeight;
  
  state.three.camera.aspect = width / height;
  state.three.camera.updateProjectionMatrix();
  state.three.renderer.setSize(width, height);
}

function animate() {
  requestAnimationFrame(animate);
  
  if (state.three.controls) {
    state.three.controls.update();
  }
  
  if (state.three.cryostat && !state.isSimulating) {
    state.three.cryostat.rotation.z += 0.0003;
  }
  
  animateDrift();

  // Render scene
  if (state.three.renderer && state.three.scene && state.three.camera) {
    state.three.renderer.render(state.three.scene, state.three.camera);
  }
}




function generateParticleHits(type) {
  const hits = [];
  const g = state.detectorGeometry;
  
  const startX = (Math.random() - 0.5) * 0.4 + 0.4; 
  const startY = (Math.random() - 0.5) * 0.4;
  const startZ = -3.5 + (Math.random() - 0.5) * 0.5;
  const vertex = new THREE.Vector3(startX, startY, startZ);
  
  state.simulationData.vertex = vertex;
  
  if (type === 'muon') {
    // 1. MUON: Long, clean, penetrating straight track.
    const length = 8.0;
    const dir = new THREE.Vector3(-0.1, -0.05, 1.0).normalize();
    const steps = 90;
    const stepSize = length / steps;
    
    for (let i = 0; i <= steps; i++) {
      const pos = vertex.clone().addScaledVector(dir, i * stepSize);
      if (isInsideTPC(pos)) {
        let charge = 0.12 + Math.random() * 0.04;
        
        if (i === 35) {
          const deltaDir = new THREE.Vector3(-0.6, 0.4, 0.2).normalize();
          for (let j = 1; j <= 8; j++) {
            const deltaPos = pos.clone().addScaledVector(deltaDir, j * 0.08);
            if (isInsideTPC(deltaPos)) {
              hits.push({ pos: deltaPos, charge: 0.15 - j * 0.015, stepIndex: i + j / 10 });
            }
          }
        }
        
        hits.push({ pos, charge, stepIndex: i });
      }
    }
  }
  else if (type === 'electron') {
    // 2. ELECTRON: Electromagnetic shower, branching tree-like structure
    generateEMShower(vertex, new THREE.Vector3(-0.08, 0.0, 1.0).normalize(), 2.8, 0, hits);
  }
  else if (type === 'photon') {
    // 3. PHOTON: Neutral conversion gap, then e+/e- shower
    const gapLength = 1.6;
    const dir = new THREE.Vector3(-0.06, 0.02, 1.0).normalize();
    const convVertex = vertex.clone().addScaledVector(dir, gapLength);
    
    const gapSteps = 15;
    for (let i = 0; i < gapSteps; i++) {
      const pos = vertex.clone().addScaledVector(dir, i * (gapLength / gapSteps));
      hits.push({ pos, charge: 0.0, stepIndex: i * 0.5, isPhotonGap: true });
    }
    
    const e1Dir = dir.clone().applyAxisAngle(new THREE.Vector3(0, 1, 0), 0.18).normalize();
    const e2Dir = dir.clone().applyAxisAngle(new THREE.Vector3(0, 1, 0), -0.18).normalize();
    
    generateEMShower(convVertex, e1Dir, 2.3, 0, hits);
    generateEMShower(convVertex, e2Dir, 2.0, 0, hits);
  }
  else if (type === 'proton') {
    // 4. PROTON: Heavy, short, high ionizing straight track ending with Bragg Peak
    const length = 1.6;
    const dir = new THREE.Vector3(-0.2, 0.5, 0.8).normalize();
    const steps = 35;
    const stepSize = length / steps;
    
    for (let i = 0; i <= steps; i++) {
      const pos = vertex.clone().addScaledVector(dir, i * stepSize);
      if (isInsideTPC(pos)) {
        const pathFraction = i / steps;
        let charge = 0.15 + 0.65 * Math.pow(pathFraction, 4.0);
        if (i === steps) charge = 0.9;
        
        hits.push({ pos, charge, stepIndex: i });
      }
    }
  }
  else if (type === 'pion') {
    // 5. PION: Hadronic interaction. Bending/scattering.
    const length1 = 2.5;
    const dir1 = new THREE.Vector3(0.0, -0.3, 1.0).normalize();
    const steps1 = 30;
    const stepSize1 = length1 / steps1;
    
    let endPos = vertex.clone();
    
    for (let i = 0; i <= steps1; i++) {
      const pos = vertex.clone().addScaledVector(dir1, i * stepSize1);
      if (isInsideTPC(pos)) {
        hits.push({ pos, charge: 0.16 + Math.random()*0.05, stepIndex: i });
        endPos = pos;
      }
    }
    
    const dir2 = new THREE.Vector3(0.4, 0.5, 0.7).normalize();
    const length2 = 1.8;
    const steps2 = 20;
    for (let i = 1; i <= steps2; i++) {
      const pos = endPos.clone().addScaledVector(dir2, i * (length2 / steps2));
      if (isInsideTPC(pos)) {
        hits.push({ pos, charge: 0.13 + Math.random()*0.03, stepIndex: steps1 + i });
      }
    }
    
    const dir3 = new THREE.Vector3(-0.5, -0.4, 0.5).normalize();
    const length3 = 1.2;
    const steps3 = 15;
    for (let i = 1; i <= steps3; i++) {
      const pos = endPos.clone().addScaledVector(dir3, i * (length3 / steps3));
      if (isInsideTPC(pos)) {
        const q = 0.25 + 0.15 * (i / steps3);
        hits.push({ pos, charge: q, stepIndex: steps1 + i });
      }
    }
  }
  
  hits.sort((a, b) => a.stepIndex - b.stepIndex);
  state.simulationData.hits = hits;
  state.simulationData.particleType = type;
  return hits;
}

function generateEMShower(start, dir, length, depth, hits) {
  if (length < 0.2 || depth > 4) return;
  
  const steps = Math.ceil(length * 15);
  const stepSize = length / steps;
  let currentPos = start.clone();
  
  for (let i = 1; i <= steps; i++) {
    currentPos = start.clone().addScaledVector(dir, i * stepSize);
    if (!isInsideTPC(currentPos)) break;
    
    const charge = (0.22 - depth * 0.04) * (1.0 - (i / steps) * 0.4) + Math.random() * 0.05;
    hits.push({ pos: currentPos.clone(), charge, stepIndex: depth * 10 + i });
  }
  
  if (depth < 4) {
    const numBranches = Math.random() > 0.4 ? 2 : 1;
    for (let b = 0; b < numBranches; b++) {
      const angle = 0.15 + Math.random() * 0.25;
      const phi = Math.random() * Math.PI * 2;
      
      const newDir = dir.clone();
      const u = new THREE.Vector3(0, 1, 0).cross(newDir).normalize();
      if (u.lengthSq() < 0.01) u.set(1, 0, 0);
      const v = newDir.clone().cross(u).normalize();
      
      newDir.addScaledVector(u, Math.sin(angle) * Math.cos(phi));
      newDir.addScaledVector(v, Math.sin(angle) * Math.sin(phi));
      newDir.normalize();
      
      generateEMShower(currentPos, newDir, length * (0.5 + Math.random() * 0.35), depth + 1, hits);
    }
  }
}

function isInsideTPC(v3) {
  const g = state.detectorGeometry;
  return (
    v3.x >= -g.widthX/2 && v3.x <= g.widthX/2 &&
    v3.y >= -g.heightY/2 && v3.y <= g.heightY/2 &&
    v3.z >= -g.lengthZ/2 && v3.z <= g.lengthZ/2
  );
}


function runSimulationEvent() {
  if (state.isSimulating) return;
  state.isSimulating = true;
  
  const btn = document.getElementById('btn-simulate');
  btn.disabled = true;
  btn.innerHTML = `<span class="glow-text-cyan">⚙</span> Simulating...`;
  
  clearSimulationScene();
  clearReadouts();
  
  const particle = state.selectedParticle;
  logMessage(`Event started. Generating interaction vertex...`);
  
  const hits = generateParticleHits(particle);
  
  logMessage(`Phase 1: Energy Deposition (Track Formation) simulated.`, true);
  animateTrackDrawing(hits, 0, () => {
    logMessage(`Phase 2: Ionization electron drift (Charge Drift) active.`, true);
    triggerChargeDrift(hits, () => {
      logMessage(`Simulation completed successfully. U, V, Y wire projections captured.`, true);
      
      btn.disabled = false;
      btn.innerHTML = `<span class="glow-text-cyan">✦</span> Simulate Event`;
      state.isSimulating = false;
    });
  });
}

function clearSimulationScene() {
  while(state.three.trackMeshGroup.children.length > 0){
    const obj = state.three.trackMeshGroup.children[0];
    state.three.trackMeshGroup.remove(obj);
  }
  
  state.three.driftParticles.forEach(p => {
    state.three.scene.remove(p.mesh);
  });
  state.three.driftParticles = [];
}


function animateTrackDrawing(hits, index, onComplete) {
  if (index >= hits.length) {
    onComplete();
    return;
  }
  
  const batchSize = Math.max(1, Math.ceil(hits.length / 30));
  const nextIndex = Math.min(hits.length, index + batchSize);
  
  for (let i = index; i < nextIndex; i++) {
    const hit = hits[i];
    
    let color = getParticleHexColor(state.selectedParticle);
    let size = 0.035;
    
    if (hit.isPhotonGap) {
      color = 0x8888ff;
      size = 0.015;
      
      const dotGeo = new THREE.SphereGeometry(size, 4, 4);
      const dotMat = new THREE.MeshBasicMaterial({ color, transparent: true, opacity: 0.5 });
      const dot = new THREE.Mesh(dotGeo, dotMat);
      dot.position.copy(hit.pos);
      state.three.trackMeshGroup.add(dot);
      continue;
    }
    
    const sphereGeo = new THREE.SphereGeometry(size + hit.charge * 0.05, 8, 8);
    const sphereMat = new THREE.MeshStandardMaterial({
      color,
      emissive: color,
      emissiveIntensity: 0.5,
      roughness: 0.2,
      metalness: 0.2
    });
    
    const sphere = new THREE.Mesh(sphereGeo, sphereMat);
    sphere.position.copy(hit.pos);
    state.three.trackMeshGroup.add(sphere);
  }
  
  setTimeout(() => {
    animateTrackDrawing(hits, nextIndex, onComplete);
  }, 35);
}


function triggerChargeDrift(hits, onComplete) {
  const g = state.detectorGeometry;
  
  const ionizingHits = hits.filter(h => !h.isPhotonGap);
  
  ionizingHits.forEach(hit => {
    const particleGeo = new THREE.SphereGeometry(0.02, 4, 4);
    const particleMat = new THREE.MeshBasicMaterial({
      color: 0xea580c,
      transparent: true,
      opacity: 0.95
    });
    
    const mesh = new THREE.Mesh(particleGeo, particleMat);
    mesh.position.copy(hit.pos);
    state.three.scene.add(mesh);
    
    state.three.driftParticles.push({
      mesh,
      initialPos: hit.pos.clone(),
      charge: hit.charge,
      arrived: false
    });
  });
  
  state.driftAnim = {
    onComplete,
    startTime: Date.now()
  };
}

function animateDrift() {
  if (!state.isSimulating || !state.driftAnim || state.three.driftParticles.length === 0) return;
  
  const g = state.detectorGeometry;
  const driftVelocity = 0.05;
  
  let allArrived = true;
  
  state.three.driftParticles.forEach(p => {
    if (p.arrived) return;
    
    allArrived = false;
    
    p.mesh.position.x -= driftVelocity;
    
    if (p.mesh.position.x <= g.anodeX) {
      p.mesh.position.x = g.anodeX;
      p.arrived = true;
      p.mesh.visible = false;
      
      const colorStr = getParticleRGBAColor(state.selectedParticle, p.charge);
      drawProjectedHit(p.initialPos, p.charge, colorStr);
    }
  });
  
  if (allArrived) {
    state.three.driftParticles.forEach(p => {
      state.three.scene.remove(p.mesh);
    });
    state.three.driftParticles = [];
    
    const cb = state.driftAnim.onComplete;
    state.driftAnim = null;
    cb();
  }
}


function getParticleHexColor(type) {
  switch (type) {
    case 'electron': return 0x0891b2; // Cyan/Teal
    case 'photon':   return 0x9333ea; // Purple
    case 'muon':     return 0x059669; // Emerald
    case 'pion':     return 0xd97706; // Amber
    case 'proton':   return 0xe11d48; // Rose
    default:         return 0x2563eb; // Blue
  }
}

function getParticleRGBAColor(type, charge) {
  const alpha = Math.max(0.4, Math.min(1.0, 0.5 + charge * 1.5));
  switch (type) {
    case 'electron': return `rgba(8, 145, 178, ${alpha})`;
    case 'photon':   return `rgba(147, 51, 234, ${alpha})`;
    case 'muon':     return `rgba(5, 150, 105, ${alpha})`;
    case 'pion':     return `rgba(217, 119, 6, ${alpha})`;
    case 'proton':   return `rgba(225, 29, 72, ${alpha})`;
    default:         return `rgba(37, 99, 235, ${alpha})`;
  }
}

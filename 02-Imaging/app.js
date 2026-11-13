/* ============================================================
   MedScope · 医学影像工作站
   纯前端：Shepp-Logan 体模 → MPR 三视图 / 窗宽窗位 / 伪彩
   光线投射体绘制与 MIP / ROI 一阶影像组学特征
   ============================================================ */
'use strict';

const $ = (s, r = document) => r.querySelector(s);
const $$ = (s, r = document) => Array.from(r.querySelectorAll(s));
const clamp = (v, a, b) => (v < a ? a : v > b ? b : v);
const fmt = (n, d = 0) => Number(n).toFixed(d);

const N = 128;                 // 体数据边长
const HU_MIN = -1000;          // LUT 覆盖的 HU 下限
const HU_SPAN = 3000;          // -1000 ~ 2000

/* ---------------- 状态 ---------------- */
const state = {
  vol: null,
  cross: { x: 64, y: 64, z: 64 },
  ww: 80,
  wl: 40,
  lut: 'gray',
  lutBuf: null,
  tool: 'cross',
  render: 'surface',
  az: 35,
  el: 18,
  zoom: 1,
  roi: null,          // { cx, cy, r } 体素坐标
  cine: null,
  mipBusy: false,
};

const CROSS_COLOR = '#22d3ee';

/* ---------------- 伪彩 ---------------- */
const LUT_STOPS = {
  gray: [[0, [0, 0, 0]], [1, [255, 255, 255]]],
  bone: [
    [0, [8, 6, 5]],
    [0.35, [86, 62, 44]],
    [0.62, [176, 132, 96]],
    [0.85, [232, 208, 178]],
    [1, [255, 250, 240]],
  ],
  hot: [
    [0, [0, 0, 0]],
    [0.32, [138, 8, 8]],
    [0.58, [220, 84, 6]],
    [0.78, [252, 184, 24]],
    [0.92, [255, 240, 140]],
    [1, [255, 255, 255]],
  ],
  rainbow: [
    [0, [10, 31, 143]],
    [0.22, [0, 150, 235]],
    [0.42, [34, 224, 107]],
    [0.62, [230, 232, 0]],
    [0.82, [255, 138, 0]],
    [1, [235, 20, 0]],
  ],
};

function lutColor(name, t) {
  const stops = LUT_STOPS[name] || LUT_STOPS.gray;
  t = clamp(t, 0, 1);
  let a = stops[0], b = stops[stops.length - 1];
  for (let i = 0; i < stops.length - 1; i++) {
    if (t >= stops[i][0] && t <= stops[i + 1][0]) { a = stops[i]; b = stops[i + 1]; break; }
  }
  const k = (t - a[0]) / (b[0] - a[0] || 1);
  return [
    a[1][0] + (b[1][0] - a[1][0]) * k,
    a[1][1] + (b[1][1] - a[1][1]) * k,
    a[1][2] + (b[1][2] - a[1][2]) * k,
  ];
}

/** 预计算 3001 级 HU → RGB 查找表 */
function buildLUT() {
  const buf = new Uint8ClampedArray((HU_SPAN + 1) * 3);
  const lutName = state.lut;
  const lo = state.wl - state.ww / 2;
  const s = 255 / state.ww;
  for (let i = 0; i <= HU_SPAN; i++) {
    const hu = HU_MIN + i;
    const t = clamp((hu - lo) * s / 255, 0, 1);
    const c = lutColor(lutName, t);
    buf[i * 3] = c[0];
    buf[i * 3 + 1] = c[1];
    buf[i * 3 + 2] = c[2];
  }
  state.lutBuf = buf;
}

/* ---------------- 体数据构建 ---------------- */
function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/* 3D Shepp–Logan 椭球参数 [a, b, c, x0, y0, z0, φ(deg), A]
   外层壳略作加厚（0.69→0.705），使三维重建的颅骨面更连续 */
const ELLIPSOIDS = [
  [0.705, 0.935, 0.915, 0, 0, 0, 0, 2.0],
  [0.6624, 0.874, 0.88, 0, -0.0184, 0, 0, -0.98],
  [0.11, 0.31, 0.22, 0.22, 0, 0, -18, -0.02],
  [0.16, 0.41, 0.28, -0.22, 0, 0, 18, -0.02],
  [0.21, 0.25, 0.41, 0, 0.35, -0.15, 0, 0.01],
  [0.046, 0.046, 0.05, 0, 0.1, 0.25, 0, 0.01],
  [0.046, 0.046, 0.05, 0, -0.1, 0.25, 0, 0.01],
  [0.046, 0.023, 0.05, -0.08, -0.605, 0, 0, 0.01],
  [0.023, 0.023, 0.02, 0, -0.606, 0, 0, 0.01],
  [0.023, 0.046, 0.02, 0.06, -0.605, 0, 0, 0.01],
];

function buildVolume() {
  const vol = new Float32Array(N * N * N);
  const cosPhi = ELLIPSOIDS.map((e) => Math.cos((e[6] * Math.PI) / 180));
  const sinPhi = ELLIPSOIDS.map((e) => Math.sin((e[6] * Math.PI) / 180));

  for (let z = 0; z < N; z++) {
    const zn = (z / (N - 1)) * 2 - 1;
    for (let y = 0; y < N; y++) {
      const yn = 1 - (y / (N - 1)) * 2;
      const rowBase = (z * N + y) * N;
      for (let x = 0; x < N; x++) {
        const xn = (x / (N - 1)) * 2 - 1;
        let v = 0;
        for (let e = 0; e < ELLIPSOIDS.length; e++) {
          const el = ELLIPSOIDS[e];
          const dy = yn - el[4];
          const dx = xn - el[3];
          const rx = dx * cosPhi[e] + dy * sinPhi[e];
          const ry = -dx * sinPhi[e] + dy * cosPhi[e];
          const rz = zn - el[5];
          const q = (rx * rx) / (el[0] * el[0]) + (ry * ry) / (el[1] * el[1]) + (rz * rz) / (el[2] * el[2]);
          if (q < 1.015) {
            // 软边界：模拟部分容积效应，避免硬阈值造成的锯齿
            const s = q <= 0.985 ? 1 : (1.015 - q) / 0.03;
            v += el[7] * s;
          }
        }
        vol[rowBase + x] = v;
      }
    }
  }

  // 模拟病灶：左额叶高密度球（影像组学 ROI 的观测目标）
  const lx = 0.34, ly = 0.42, lz = 0.02, lr = 0.135;
  for (let z = 0; z < N; z++) {
    const zn = (z / (N - 1)) * 2 - 1;
    const dz = zn - lz;
    if (Math.abs(dz) > lr) continue;
    for (let y = 0; y < N; y++) {
      const yn = 1 - (y / (N - 1)) * 2;
      const dy = yn - ly;
      if (Math.abs(dy) > lr) continue;
      const rowBase = (z * N + y) * N;
      for (let x = 0; x < N; x++) {
        const xn = (x / (N - 1)) * 2 - 1;
        const dx = xn - lx;
        const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (d <= lr) vol[rowBase + x] += 0.075 * (1 - (d / lr) * 0.35);
      }
    }
  }

  // 成像噪声（±3 HU 量级）
  const rnd = mulberry32(20261113);
  for (let i = 0; i < vol.length; i++) {
    if (vol[i] > 0.2) vol[i] += (rnd() - 0.5) * 0.006;
  }
  state.vol = vol;
  state.vol3d = blurVolume(vol);
  state.volMip = stripBone(vol);
}

/** 三维重建专用的轻度平滑体数据（可分离 3 点盒式滤波） */
function blurVolume(src) {
  const n3 = N * N * N;
  let a = Float32Array.from(src);
  let b = new Float32Array(n3);
  for (let axis = 0; axis < 3; axis++) {
    for (let z = 0; z < N; z++) {
      for (let y = 0; y < N; y++) {
        const base = (z * N + y) * N;
        for (let x = 0; x < N; x++) {
          let i0, i1, i2;
          if (axis === 0) {
            const zz = z, yy = y, xx = x;
            i0 = (zz * N + yy) * N + Math.max(0, xx - 1);
            i1 = (zz * N + yy) * N + xx;
            i2 = (zz * N + yy) * N + Math.min(N - 1, xx + 1);
          } else if (axis === 1) {
            i0 = (z * N + Math.max(0, y - 1)) * N + x;
            i1 = (z * N + y) * N + x;
            i2 = (z * N + Math.min(N - 1, y + 1)) * N + x;
          } else {
            i0 = (Math.max(0, z - 1) * N + y) * N + x;
            i1 = (z * N + y) * N + x;
            i2 = (Math.min(N - 1, z + 1) * N + y) * N + x;
          }
          b[base + x] = (a[i0] + a[i1] + a[i2]) / 3;
        }
      }
    }
    const t = a; a = b; b = t;
  }
  return a;
}

/* 体素值 → HU */
const huAt = (idx) => state.vol[idx] * 1000 - 980;

/* ---------------- 画布工具 ---------------- */
function fitCanvas(cv) {
  // 始终以元素实际显示尺寸为准（CSS 100%），避免读到过期尺寸导致画面溢出
  const rect = cv.getBoundingClientRect();
  const w = Math.max(1, Math.round(rect.width));
  const h = Math.max(1, Math.round(rect.height));
  const dpr = Math.min(2, window.devicePixelRatio || 1);
  cv.width = Math.round(w * dpr);
  cv.height = Math.round(h * dpr);
  const ctx = cv.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  return { ctx, w, h, dpr };
}

const sliceCanvas = document.createElement('canvas');
sliceCanvas.width = N;
sliceCanvas.height = N;
const sliceCtx = sliceCanvas.getContext('2d');
const sliceImage = sliceCtx.createImageData(N, N);

/** 生成指定视图的二维切片像素（写入 sliceCanvas） */
function fillSlice(view) {
  const data = sliceImage.data;
  const buf = state.lutBuf;
  const vol = state.vol;
  const { x: cx, y: cy, z: cz } = state.cross;

  for (let j = 0; j < N; j++) {
    for (let i = 0; i < N; i++) {
      let idx;
      if (view === 'axial') {
        idx = (cz * N + j) * N + i;
      } else if (view === 'coronal') {
        const z = N - 1 - j;
        idx = (z * N + cy) * N + i;
      } else {
        const z = N - 1 - j;
        idx = (z * N + i) * N + cx;   // i 对应 y
      }
      const hu = vol[idx] * 1000 - 980;
      let k = hu + 1000;
      k = k < 0 ? 0 : k > HU_SPAN ? HU_SPAN : k | 0;
      const o = (j * N + i) * 4;
      const p = k * 3;
      data[o] = buf[p];
      data[o + 1] = buf[p + 1];
      data[o + 2] = buf[p + 2];
      data[o + 3] = 255;
    }
  }
  sliceCtx.putImageData(sliceImage, 0, 0);
}

const ORIENT = {
  axial: { t: 'A', b: 'P', l: 'R', r: 'L' },
  coronal: { t: 'S', b: 'I', l: 'R', r: 'L' },
  sagittal: { t: 'S', b: 'I', l: 'A', r: 'P' },
};

function drawView(view) {
  const cv = $('#cv' + view[0].toUpperCase() + view.slice(1));
  const box = cv.parentElement;
  const w = box.clientWidth, h = box.clientHeight;
  if (!w || !h) return;
  const { ctx } = fitCanvas(cv, w, h);
  ctx.fillStyle = '#000';
  ctx.fillRect(0, 0, w, h);

  const size = Math.max(40, Math.min(w, h) - 26);
  const ox = (w - size) / 2, oy = (h - size) / 2;
  ctx.imageSmoothingEnabled = true;
  ctx.imageSmoothingQuality = 'high';

  fillSlice(view);
  ctx.drawImage(sliceCanvas, ox, oy, size, size);

  // 边框
  ctx.strokeStyle = 'rgba(120,150,175,.16)';
  ctx.lineWidth = 1;
  ctx.strokeRect(ox - 0.5, oy - 0.5, size + 1, size + 1);

  const px = (v) => ox + (v / (N - 1)) * size;
  const py = (v) => oy + (v / (N - 1)) * size;

  // 十字定位线
  let lx, ly;
  if (view === 'axial') { lx = px(state.cross.x); ly = py(state.cross.y); }
  else if (view === 'coronal') { lx = px(state.cross.x); ly = py(N - 1 - state.cross.z); }
  else { lx = px(state.cross.y); ly = py(N - 1 - state.cross.z); }

  ctx.save();
  ctx.strokeStyle = 'rgba(34,211,238,.55)';
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(lx, oy); ctx.lineTo(lx, oy + size);
  ctx.moveTo(ox, ly); ctx.lineTo(ox + size, ly);
  ctx.stroke();
  ctx.strokeStyle = CROSS_COLOR;
  ctx.lineWidth = 1.4;
  ctx.beginPath();
  ctx.arc(lx, ly, 5, 0, Math.PI * 2);
  ctx.stroke();
  ctx.restore();

  // 方位标记
  const o = ORIENT[view];
  ctx.fillStyle = 'rgba(140,170,195,.72)';
  ctx.font = '600 10px "Cascadia Mono", Consolas, monospace';
  ctx.textAlign = 'center';
  ctx.fillText(o.t, ox + size / 2, oy + 12);
  ctx.fillText(o.b, ox + size / 2, oy + size - 5);
  ctx.textAlign = 'left';
  ctx.fillText(o.l, ox + 5, oy + size / 2 + 3.5);
  ctx.textAlign = 'right';
  ctx.fillText(o.r, ox + size - 5, oy + size / 2 + 3.5);

  // 比例尺 20 mm
  const barLen = (20 / 128) * size;
  ctx.fillStyle = 'rgba(160,190,210,.7)';
  ctx.fillRect(ox + 10, oy + size - 16, barLen, 2);
  ctx.textAlign = 'left';
  ctx.font = '9px "Cascadia Mono", Consolas, monospace';
  ctx.fillText('20 mm', ox + 10, oy + size - 21);

  // ROI 圆
  if (view === 'axial' && state.roi) {
    const r = state.roi;
    ctx.save();
    ctx.setLineDash([3, 3]);
    ctx.strokeStyle = '#fbbf24';
    ctx.lineWidth = 1.3;
    ctx.beginPath();
    ctx.arc(px(r.cx), py(r.cy), (r.r / (N - 1)) * size, 0, Math.PI * 2);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = 'rgba(251,191,36,.9)';
    ctx.font = '9.5px "Cascadia Mono", Consolas, monospace';
    ctx.textAlign = 'left';
    ctx.fillText(`ROI r=${r.r.toFixed(1)}mm`, px(r.cx) + 6, py(r.cy) - 6);
    ctx.restore();
  }

  // 视图信息
  if (view === 'axial') $('#posAxial').textContent = `z = ${state.cross.z + 1}`;
  if (view === 'coronal') $('#posCoronal').textContent = `y = ${state.cross.y + 1}`;
  if (view === 'sagittal') $('#posSagittal').textContent = `x = ${state.cross.x + 1}`;
}

/* ---------------- 3D 光线投射 ---------------- */
/* 体绘制传递函数：软组织 / 病灶高密度带 / 颅骨 三段 */
const DVR = {
  soft: 0.0012, softFrom: 0.6, softTo: 1.055, softColor: [0.74, 0.84, 0.94],
  focus: 0.03, focusColor: [1.0, 0.72, 0.34],
  bone: 0.13, boneFrom: 1.55, boneColor: [1, 0.99, 0.95],
};

function sampleVol(x, y, z) {
  // x,y,z 归一化 [-1,1] → 体素（三维重建使用平滑后的体数据）
  const ix = ((x + 1) * 0.5 * (N - 1)) | 0;
  const iy = ((1 - y) * 0.5 * (N - 1)) | 0;
  const iz = ((z + 1) * 0.5 * (N - 1)) | 0;
  if (ix < 0 || ix >= N || iy < 0 || iy >= N || iz < 0 || iz >= N) return 0;
  const v = state.vol3d || state.vol;
  return v[(iz * N + iy) * N + ix];
}

let mipCanvas = document.createElement('canvas');
let mipCtx = mipCanvas.getContext('2d');

const RENDER_LABELS = { surface: '表面重建', dvr: '体绘制', mip: 'MIP' };
const renderLabel = () => RENDER_LABELS[state.render] || '表面重建';

/** 去骨：对高密度区做掩膜并膨胀，抹掉颅骨及其部分容积体素，供 MIP 使用 */
function stripBone(src) {
  const n3 = src.length;
  const mask = new Uint8Array(n3);
  for (let i = 0; i < n3; i++) if (src[i] > 1.45) mask[i] = 1;
  const tmp = new Uint8Array(n3);
  for (let pass = 0; pass < 2; pass++) {
    tmp.set(mask);
    for (let z = 1; z < N - 1; z++) {
      for (let y = 1; y < N - 1; y++) {
        const base = (z * N + y) * N;
        for (let x = 1; x < N - 1; x++) {
          const i = base + x;
          if (mask[i]) continue;
          if (tmp[i - 1] || tmp[i + 1] || tmp[i - N] || tmp[i + N] || tmp[i - N * N] || tmp[i + N * N]) mask[i] = 1;
        }
      }
    }
  }
  const out = new Float32Array(n3);
  for (let i = 0; i < n3; i++) out[i] = mask[i] ? 0 : src[i];
  return out;
}

/** 原始（未平滑）体数据采样：用于需要锐利边界的去骨 MIP */
function sampleRaw(x, y, z) {
  const ix = ((x + 1) * 0.5 * (N - 1)) | 0;
  const iy = ((1 - y) * 0.5 * (N - 1)) | 0;
  const iz = ((z + 1) * 0.5 * (N - 1)) | 0;
  if (ix < 0 || ix >= N || iy < 0 || iy >= N || iz < 0 || iz >= N) return 0;
  const v = state.volMip || state.vol;
  return v[(iz * N + iy) * N + ix];
}

/** 三线性插值采样（用于等值面细化与梯度计算，消除阶梯） */
function sampleTri(x, y, z) {
  const v = state.vol3d || state.vol;
  let fx = (x + 1) * 0.5 * (N - 1);
  let fy = (1 - y) * 0.5 * (N - 1);
  let fz = (z + 1) * 0.5 * (N - 1);
  if (fx < 0 || fx > N - 1 || fy < 0 || fy > N - 1 || fz < 0 || fz > N - 1) return 0;
  const ix = Math.min(N - 2, fx | 0), iy = Math.min(N - 2, fy | 0), iz = Math.min(N - 2, fz | 0);
  const tx = fx - ix, ty = fy - iy, tz = fz - iz;
  const i000 = (iz * N + iy) * N + ix;
  const i100 = i000 + 1;
  const i010 = i000 + N;
  const i110 = i010 + 1;
  const i001 = i000 + N * N;
  const i101 = i001 + 1;
  const i011 = i001 + N;
  const i111 = i011 + 1;
  const c00 = v[i000] * (1 - tx) + v[i100] * tx;
  const c10 = v[i010] * (1 - tx) + v[i110] * tx;
  const c01 = v[i001] * (1 - tx) + v[i101] * tx;
  const c11 = v[i011] * (1 - tx) + v[i111] * tx;
  const c0 = c00 * (1 - ty) + c10 * ty;
  const c1 = c01 * (1 - ty) + c11 * ty;
  return c0 * (1 - tz) + c1 * tz;
}

function render3D(quality) {
  if (!state.vol) return;
  const S = quality ? 232 : 116;
  const steps = quality ? 190 : 110;
  mipCanvas.width = S;
  mipCanvas.height = S;
  const img = mipCtx.createImageData(S, S);
  const data = img.data;
  const lut = state.lutBuf;

  const az = (state.az * Math.PI) / 180;
  const el = (state.el * Math.PI) / 180;
  const ca = Math.cos(az), sa = Math.sin(az), ce = Math.cos(el), se = Math.sin(el);

  // 相机方向（由相机指向原点）
  const dx = -ce * ca, dy = -ce * sa, dz = -se;
  // 右向量 = normalize(cross(d, up))，up = (0,0,1)
  let rx = dy * 1 - dz * 0, ry = dz * 0 - dx * 1, rz = 0;
  const rl = Math.hypot(rx, ry, rz) || 1;
  rx /= rl; ry /= rl; rz /= rl;
  // 上向量 = cross(right, d)
  const ux = ry * dz - rz * dy, uy = rz * dx - rx * dz, uz = rx * dy - ry * dx;

  const span = 1.28 * state.zoom;
  const step = (2.6 * 1.28) / steps;
  const startD = -1.35 * 1.28;
  const mode = state.render;
  const isMip = mode === 'mip';
  const isSurface = mode === 'surface';
  const THR = 1.55;

  // 头灯方向（相机左上方），保证任意角度都有良好照明
  let lx = -dx + rx * 0.42 + ux * 0.56;
  let ly = -dy + ry * 0.42 + uy * 0.56;
  let lz = -dz + rz * 0.42 + uz * 0.56;
  const ll = Math.hypot(lx, ly, lz) || 1;
  lx /= ll; ly /= ll; lz /= ll;

  for (let py = 0; py < S; py++) {
    const v = (1 - (py + 0.5) / S * 2) * span;
    for (let px = 0; px < S; px++) {
      const u = ((px + 0.5) / S * 2 - 1) * span;
      let ox = rx * u + ux * v + dx * startD;
      let oy = ry * u + uy * v + dy * startD;
      let oz = rz * u + uz * v + dz * startD;
      const sx = dx * step, sy = dy * step, sz = dz * step;

      let R = 0, G = 0, B = 0;

      if (isSurface) {
        let hx = 0, hy = 0, hz = 0, hit = false;
        for (let s = 0; s < steps; s++) {
          const px0 = ox, py0 = oy, pz0 = oz;
          ox += sx; oy += sy; oz += sz;
          if (ox > 1.02 || ox < -1.02 || oy > 1.02 || oy < -1.02 || oz > 1.02 || oz < -1.02) continue;
          if (sampleVol(ox, oy, oz) > THR) {
            // 三线性插值 + 二分细化交点，消除阶梯
            let t0 = 0, t1 = 1;
            for (let b = 0; b < 6; b++) {
              const tm = (t0 + t1) * 0.5;
              const xm = px0 + (ox - px0) * tm;
              const ym = py0 + (oy - py0) * tm;
              const zm = pz0 + (oz - pz0) * tm;
              if (sampleTri(xm, ym, zm) > THR) t1 = tm; else t0 = tm;
            }
            hx = px0 + (ox - px0) * t1;
            hy = py0 + (oy - py0) * t1;
            hz = pz0 + (oz - pz0) * t1;
            hit = true;
            break;
          }
        }
        if (hit) {
          const h = 0.028;
          const gx = sampleTri(hx + h, hy, hz) - sampleTri(hx - h, hy, hz);
          const gy = sampleTri(hx, hy + h, hz) - sampleTri(hx, hy - h, hz);
          const gz = sampleTri(hx, hy, hz + h) - sampleTri(hx, hy, hz - h);
          const gl = Math.hypot(gx, gy, gz) || 1;
          const nx = -gx / gl, ny = -gy / gl, nz = -gz / gl;
          const diff = Math.max(0, nx * lx + ny * ly + nz * lz);
          const rim = Math.max(0, -(nx * dx + ny * dy + nz * dz));
          const shade = 0.13 + 0.8 * diff + 0.2 * Math.pow(rim, 2.4);
          R = clamp(226 * shade + 40 * Math.pow(diff, 14), 0, 255);
          G = clamp(231 * shade + 40 * Math.pow(diff, 14), 0, 255);
          B = clamp(238 * shade + 40 * Math.pow(diff, 14), 0, 255);
          R += 20 * Math.pow(rim, 5);
          G += 34 * Math.pow(rim, 5);
          B += 46 * Math.pow(rim, 5);
        }
      } else {
        let acc = 0, cr = 0, cg = 0, cb = 0, maxV = 0;
        for (let s = 0; s < steps; s++) {
          ox += sx; oy += sy; oz += sz;
          if (ox > 1.05 || ox < -1.05 || oy > 1.05 || oy < -1.05 || oz > 1.05 || oz < -1.05) continue;
          const val = isMip ? sampleRaw(ox, oy, oz) : sampleVol(ox, oy, oz);
          if (isMip) {
            // 去骨 MIP：体数据已抹除颅骨，直接取最大密度
            if (val > maxV) maxV = val;
          } else {
            let a = 0, cr0 = 1, cg0 = 1, cb0 = 1;
            if (val > DVR.boneFrom) {
              a = DVR.bone;
              cr0 = DVR.boneColor[0]; cg0 = DVR.boneColor[1]; cb0 = DVR.boneColor[2];
            } else if (val > DVR.softTo) {
              a = DVR.focus;
              cr0 = DVR.focusColor[0]; cg0 = DVR.focusColor[1]; cb0 = DVR.focusColor[2];
            } else if (val > DVR.softFrom) {
              a = DVR.soft;
              cr0 = DVR.softColor[0]; cg0 = DVR.softColor[1]; cb0 = DVR.softColor[2];
            }
            if (a > 0) {
              const wgt = a * (1 - acc);
              cr += cr0 * wgt;
              cg += cg0 * wgt;
              cb += cb0 * wgt;
              acc += wgt;
              if (acc > 0.985) break;
            }
          }
        }

        if (isMip) {
          const hu = maxV * 1000 - 980;
          let k = hu + 1000;
          k = k < 0 ? 0 : k > HU_SPAN ? HU_SPAN : k | 0;
          R = lut[k * 3]; G = lut[k * 3 + 1]; B = lut[k * 3 + 2];
        } else {
          R = clamp(cr, 0, 1) * 255;
          G = clamp(cg, 0, 1) * 255;
          B = clamp(cb, 0, 1) * 255;
        }
      }
      const o = (py * S + px) * 4;
      data[o] = R; data[o + 1] = G; data[o + 2] = B; data[o + 3] = 255;
    }
  }
  mipCtx.putImageData(img, 0, 0);

  // 绘制到视口
  const cv = $('#cvMip');
  const box = cv.parentElement;
  const w = box.clientWidth, h = box.clientHeight;
  if (!w || !h) return;
  const { ctx } = fitCanvas(cv, w, h);
  ctx.fillStyle = '#000';
  ctx.fillRect(0, 0, w, h);
  const size = Math.min(w, h) - 20;
  ctx.imageSmoothingEnabled = true;
  ctx.imageSmoothingQuality = 'high';
  ctx.drawImage(mipCanvas, (w - size) / 2, (h - size) / 2, size, size);

  // 轴指示
  ctx.save();
  ctx.translate(w - 46, h - 46);
  const axes = [
    { n: 'A', v: [0, 1, 0], c: '#22d3ee' },
    { n: 'L', v: [1, 0, 0], c: '#fbbf24' },
    { n: 'S', v: [0, 0, 1], c: '#a78bfa' },
  ];
  for (const ax of axes) {
    const X = ax.v[0] * ca + ax.v[1] * sa;
    const Y = (-ax.v[0] * sa + ax.v[1] * ca) * se + ax.v[2] * ce;
    ctx.strokeStyle = ax.c;
    ctx.globalAlpha = 0.8;
    ctx.lineWidth = 1.2;
    ctx.beginPath();
    ctx.moveTo(0, 0);
    ctx.lineTo(X * 20, -Y * 20);
    ctx.stroke();
    ctx.globalAlpha = 1;
    ctx.fillStyle = ax.c;
    ctx.font = '9px "Cascadia Mono", Consolas, monospace';
    ctx.fillText(ax.n, X * 26 - 3, -Y * 26 + 3);
  }
  ctx.restore();

  $('#mipAngles').textContent = `az ${state.az.toFixed(0)}° / el ${state.el.toFixed(0)}°`;
  $('#mipInfo').textContent = `${renderLabel()} · ${steps} steps`;
}

/* ---------------- ROI 影像组学 ---------------- */
function roiStats() {
  const { cx, cy, r } = state.roi || {};
  if (cx == null) return null;
  const z = state.cross.z;
  const vals = [];
  const x0 = Math.max(0, Math.floor(cx - r)), x1 = Math.min(N - 1, Math.ceil(cx + r));
  const y0 = Math.max(0, Math.floor(cy - r)), y1 = Math.min(N - 1, Math.ceil(cy + r));
  for (let y = y0; y <= y1; y++) {
    for (let x = x0; x <= x1; x++) {
      const dx = x - cx, dy = y - cy;
      if (dx * dx + dy * dy > r * r) continue;
      vals.push(huAt((z * N + y) * N + x));
    }
  }
  if (vals.length < 4) return null;
  const n = vals.length;
  let sum = 0;
  for (const v of vals) sum += v;
  const mean = sum / n;
  let m2 = 0, m3 = 0, m4 = 0;
  for (const v of vals) {
    const d = v - mean;
    m2 += d * d; m3 += d * d * d; m4 += d * d * d * d;
  }
  const variance = m2 / n;
  const std = Math.sqrt(variance);
  const skew = std ? (m3 / n) / (std * std * std) : 0;
  const kurt = variance ? (m4 / n) / (variance * variance) - 3 : 0;
  const sorted = Float64Array.from(vals).sort();
  const median = n % 2 ? sorted[(n - 1) / 2] : (sorted[n / 2 - 1] + sorted[n / 2]) / 2;
  const p10 = sorted[Math.floor(n * 0.1)];
  const p90 = sorted[Math.floor(n * 0.9)];

  // 32 bin 直方图 + 香农熵
  const lo = sorted[0], hi = sorted[n - 1];
  const bins = 32;
  const hist = new Float64Array(bins);
  const bw = (hi - lo) / bins || 1;
  for (const v of vals) hist[Math.min(bins - 1, Math.floor((v - lo) / bw))]++;
  let entropy = 0;
  for (let i = 0; i < bins; i++) {
    const p = hist[i] / n;
    if (p > 0) entropy -= p * Math.log2(p);
  }

  return {
    n, mean, std, min: sorted[0], max: sorted[n - 1], median, p10, p90,
    skew, kurt, entropy, lo, hi, hist, bins,
    area: (n * 1.0).toFixed(0),
  };
}

function renderRoi() {
  const st = roiStats();
  if (!st) {
    $('#roiEmpty').classList.remove('hidden');
    $('#roiPanel').classList.add('hidden');
    return;
  }
  $('#roiEmpty').classList.add('hidden');
  $('#roiPanel').classList.remove('hidden');
  $('#roiCenter').textContent = `(${state.roi.cx.toFixed(0)}, ${state.roi.cy.toFixed(0)})`;
  $('#roiRadius').textContent = `${state.roi.r.toFixed(1)} mm`;

  // 直方图
  const cv = $('#cvHist');
  const { ctx, w, h } = fitCanvas(cv, cv.parentElement.clientWidth, 92);
  ctx.clearRect(0, 0, w, h);
  const pad = 4;
  const bw = (w - pad * 2) / st.bins;
  const maxBin = Math.max(...st.hist);
  for (let i = 0; i < st.bins; i++) {
    const bh = (st.hist[i] / maxBin) * (h - pad * 2);
    ctx.fillStyle = i % 2 ? 'rgba(34,211,238,.5)' : 'rgba(34,211,238,.72)';
    ctx.fillRect(pad + i * bw, h - pad - bh, Math.max(1, bw - 1), bh);
  }
  ctx.strokeStyle = 'rgba(120,150,175,.25)';
  ctx.beginPath();
  ctx.moveTo(pad, h - pad);
  ctx.lineTo(w - pad, h - pad);
  ctx.stroke();
  ctx.fillStyle = 'rgba(140,170,195,.65)';
  ctx.font = '8.5px "Cascadia Mono", Consolas, monospace';
  ctx.textAlign = 'left';
  ctx.fillText(`${st.lo.toFixed(0)} HU`, pad, 9);
  ctx.textAlign = 'right';
  ctx.fillText(`${st.hi.toFixed(0)} HU`, w - pad, 9);

  const rows = [
    ['面积', `${st.n} mm²`],
    ['均值 Mean', `${st.mean.toFixed(1)} HU`],
    ['标准差 Std', `${st.std.toFixed(1)} HU`],
    ['中位数', `${st.median.toFixed(1)} HU`],
    ['P10 / P90', `${st.p10.toFixed(0)} / ${st.p90.toFixed(0)}`],
    ['偏度 Skew', st.skew.toFixed(3)],
    ['峰度 Kurt', st.kurt.toFixed(3)],
    ['均匀性 Entropy', st.entropy.toFixed(3)],
    ['最小值 / 最大值', `${st.min.toFixed(0)} / ${st.max.toFixed(0)}`],
  ];
  $('#featList').innerHTML = rows
    .map(([k, v]) => `<div class="feat"><span>${k}</span><b>${v}</b></div>`)
    .join('');
}

function demoRoi() {
  const cx = Math.round(((0.34 + 1) / 2) * (N - 1));
  const cy = Math.round(((1 - 0.42) / 2) * (N - 1));
  const cz = Math.round(((0.02 + 1) / 2) * (N - 1));
  state.cross.z = cz;
  state.roi = { cx, cy, r: 13 };
  state.tool = 'roi';
  $('#toolRoi').classList.add('on');
  $('#toolCross').classList.remove('on');
  syncSliceUI();
  updateAll();
}

/* ---------------- 胶片条 ---------------- */
function buildFilm() {
  const film = $('#film');
  film.innerHTML = '';
  const zs = [6, 22, 38, 54, 70, 86, 102, 118];
  for (const z of zs) {
    const btn = document.createElement('button');
    const cv = document.createElement('canvas');
    cv.width = 64; cv.height = 64;
    const c = cv.getContext('2d');
    const tmp = c.createImageData(64, 64);
    const buf = state.lutBuf;
    for (let j = 0; j < 64; j++) {
      const sy = Math.min(N - 1, Math.floor((j / 64) * N));
      for (let i = 0; i < 64; i++) {
        const sx = Math.min(N - 1, Math.floor((i / 64) * N));
        const hu = huAt((z * N + sy) * N + sx);
        let k = clamp(hu + 1000, 0, HU_SPAN) | 0;
        const o = (j * 64 + i) * 4;
        tmp.data[o] = buf[k * 3];
        tmp.data[o + 1] = buf[k * 3 + 1];
        tmp.data[o + 2] = buf[k * 3 + 2];
        tmp.data[o + 3] = 255;
      }
    }
    c.putImageData(tmp, 0, 0);
    const label = document.createElement('span');
    label.textContent = z + 1;
    btn.appendChild(cv);
    btn.appendChild(label);
    btn.dataset.z = z;
    btn.addEventListener('click', () => {
      state.cross.z = z;
      syncSliceUI();
      updateAll();
    });
    film.appendChild(btn);
  }
  markFilm();
}

function markFilm() {
  const cur = state.cross.z;
  $$('#film button').forEach((b) => {
    const z = Number(b.dataset.z);
    b.classList.toggle('on', Math.abs(z - cur) <= 8);
  });
}

/* ---------------- 更新流程 ---------------- */
let rafId = null;
let mipTimer = null;

function updateAll() {
  if (rafId) return;
  rafId = requestAnimationFrame(() => {
    rafId = null;
    drawView('axial');
    drawView('coronal');
    drawView('sagittal');
    renderRoi();
    markFilm();
    $('#probeCross').textContent = `${state.cross.x + 1} / ${state.cross.y + 1} / ${state.cross.z + 1}`;
  });
}

function update3DPreview() {
  clearTimeout(mipTimer);
  render3D(false);
}

/* ---------------- 交互 ---------------- */
function syncSliceUI() {
  $('#sliceRange').value = state.cross.z;
  $('#sliceNum').textContent = state.cross.z + 1;
}

function viewCoords(view, ev) {
  const cv = $('#cv' + view[0].toUpperCase() + view.slice(1));
  const rect = cv.getBoundingClientRect();
  const w = rect.width, h = rect.height;
  const size = Math.max(40, Math.min(w, h) - 26);
  const ox = (w - size) / 2, oy = (h - size) / 2;
  const x = ev.clientX - rect.left, y = ev.clientY - rect.top;
  const u = (x - ox) / size, v = (y - oy) / size;
  if (u < 0 || u > 1 || v < 0 || v > 1) return null;
  const i = clamp(Math.round(u * (N - 1)), 0, N - 1);
  const j = clamp(Math.round(v * (N - 1)), 0, N - 1);
  return { i, j, u, v };
}

function pickVoxel(view, ev) {
  const p = viewCoords(view, ev);
  if (!p) return null;
  const c = { ...state.cross };
  if (view === 'axial') { c.x = p.i; c.y = p.j; }
  else if (view === 'coronal') { c.x = p.i; c.z = N - 1 - p.j; }
  else { c.y = p.i; c.z = N - 1 - p.j; }
  return c;
}

let drag = null;

function bindViewports() {
  // 视口尺寸变化时重绘，避免布局抖动导致画布与容器尺寸不一致
  const sizeTimers = {};
  if (window.ResizeObserver) {
    const ro = new ResizeObserver((entries) => {
      for (const en of entries) {
        const view = en.target.dataset.view;
        if (view === 'mip') {
          clearTimeout(sizeTimers.mip);
          sizeTimers.mip = setTimeout(() => render3D(true), 180);
        } else {
          clearTimeout(sizeTimers[view]);
          sizeTimers[view] = setTimeout(() => drawView(view), 60);
        }
      }
    });
    $$('.vp').forEach((vp) => ro.observe(vp));
  }

  ['axial', 'coronal', 'sagittal'].forEach((view) => {
    const vp = $(`.vp[data-view="${view}"]`);
    const cv = $('#cv' + view[0].toUpperCase() + view.slice(1));
    const huLabel = $('#hu' + view[0].toUpperCase() + view.slice(1));

    cv.addEventListener('mousemove', (ev) => {
      const c = pickVoxel(view, ev);
      if (!c) { huLabel.textContent = '—'; $('#tooltip').classList.remove('show'); return; }
      const idx = (c.z * N + c.y) * N + c.x;
      const hu = huAt(idx);
      huLabel.textContent = `${hu.toFixed(0)} HU`;
      $('#probeHU').textContent = `${hu.toFixed(0)} HU`;
      $('#probeXYZ').textContent = `${c.x + 1} / ${c.y + 1} / ${c.z + 1}`;

      if (drag && drag.type === 'cross' && drag.view === view) {
        state.cross.x = c.x; state.cross.y = c.y; state.cross.z = c.z;
        syncSliceUI();
        updateAll();
      } else if (drag && drag.type === 'roi' && drag.view === 'axial') {
        const dx = c.x - drag.cx, dy = c.y - drag.cy;
        state.roi = { cx: drag.cx, cy: drag.cy, r: Math.max(2, Math.hypot(dx, dy)) };
        updateAll();
      }
    });

    cv.addEventListener('mouseleave', () => {
      huLabel.textContent = '—';
      $('#tooltip').classList.remove('show');
    });

    cv.addEventListener('mousedown', (ev) => {
      const c = pickVoxel(view, ev);
      if (!c) return;
      ev.preventDefault();
      if (state.tool === 'roi' && view === 'axial') {
        drag = { type: 'roi', view, cx: c.x, cy: c.y };
        state.roi = { cx: c.x, cy: c.y, r: 2 };
      } else {
        drag = { type: 'cross', view };
        state.cross.x = c.x; state.cross.y = c.y; state.cross.z = c.z;
        syncSliceUI();
      }
      cv.parentElement.classList.add('active');
      updateAll();
    });

    cv.addEventListener('wheel', (ev) => {
      ev.preventDefault();
      const dir = ev.deltaY > 0 ? 1 : -1;
      state.cross.z = clamp(state.cross.z + dir, 0, N - 1);
      syncSliceUI();
      updateAll();
    }, { passive: false });

    // 双击：把该点设为十字中心
    cv.addEventListener('dblclick', (ev) => {
      const c = pickVoxel(view, ev);
      if (c) { state.cross.x = c.x; state.cross.y = c.y; state.cross.z = c.z; syncSliceUI(); updateAll(); }
    });
  });

  window.addEventListener('mouseup', () => {
    if (drag) {
      drag = null;
      $$('.vp').forEach((v) => v.classList.remove('active'));
    }
  });

  // 3D 视图拖拽
  const cv3 = $('#cvMip');
  let rot = null;
  cv3.addEventListener('mousedown', (ev) => {
    ev.preventDefault();
    rot = { x: ev.clientX, y: ev.clientY, az: state.az, el: state.el };
    cv3.parentElement.classList.add('active');
  });
  window.addEventListener('mousemove', (ev) => {
    if (!rot) return;
    state.az = rot.az + (ev.clientX - rot.x) * 0.5;
    state.el = clamp(rot.el + (ev.clientY - rot.y) * 0.4, -85, 85);
    update3DPreview();
  });
  window.addEventListener('mouseup', () => {
    if (rot) {
      rot = null;
      cv3.parentElement.classList.remove('active');
      render3D(true);
    }
  });
  cv3.addEventListener('wheel', (ev) => {
    ev.preventDefault();
    state.zoom = clamp(state.zoom * (ev.deltaY > 0 ? 0.92 : 1.08), 0.55, 2.4);
    update3DPreview();
    clearTimeout(mipTimer);
    mipTimer = setTimeout(() => render3D(true), 120);
  }, { passive: false });
}

/* ---------------- 控件绑定 ---------------- */
function setWindow(ww, wl, preset) {
  state.ww = ww;
  state.wl = wl;
  $('#rngWW').value = ww;
  $('#rngWL').value = wl;
  $('#outWW').textContent = ww;
  $('#outWL').textContent = wl;
  buildLUT();
  if (preset) {
    $$('.presets [data-preset]').forEach((b) => b.classList.toggle('on', b.dataset.preset === preset));
  } else {
    $$('.presets [data-preset]').forEach((b) => b.classList.remove('on'));
  }
  updateAll();
}

function bindControls() {
  $$('.presets [data-preset]').forEach((b) => {
    b.addEventListener('click', () => {
      const p = b.dataset.preset;
      const presets = {
        brain: [80, 40],
        bone: [2000, 500],
        soft: [350, 40],
        full: [3000, 300],
      };
      const [ww, wl] = presets[p];
      setWindow(ww, wl, p);
    });
  });

  $('#rngWW').addEventListener('input', () => setWindow(Number($('#rngWW').value), state.wl, null));
  $('#rngWL').addEventListener('input', () => setWindow(state.ww, Number($('#rngWL').value), null));

  $$('#luts button').forEach((b) => {
    b.addEventListener('click', () => {
      state.lut = b.dataset.lut;
      $$('#luts button').forEach((x) => x.classList.toggle('on', x === b));
      buildLUT();
      buildFilm();
      updateAll();
    });
  });

  $('#toolCross').addEventListener('click', () => {
    state.tool = 'cross';
    $('#toolCross').classList.add('on');
    $('#toolRoi').classList.remove('on');
  });
  $('#toolRoi').addEventListener('click', () => {
    state.tool = 'roi';
    $('#toolRoi').classList.add('on');
    $('#toolCross').classList.remove('on');
  });

  $$('.presets [data-render]').forEach((b) => {
    b.addEventListener('click', () => {
      state.render = b.dataset.render;
      $$('.presets [data-render]').forEach((x) => x.classList.toggle('on', x === b));
      $('#renderMode').textContent = renderLabel();
      render3D(true);
    });
  });

  $('#btnReset').addEventListener('click', () => {
    state.cross = { x: 64, y: 64, z: 64 };
    state.az = 35; state.el = 18; state.zoom = 1;
    state.roi = null;
    setWindow(80, 40, 'brain');
    syncSliceUI();
    updateAll();
  });

  $('#btnCine').addEventListener('click', (e) => {
    if (state.cine) {
      clearInterval(state.cine);
      state.cine = null;
      e.currentTarget.classList.remove('on');
    } else {
      e.currentTarget.classList.add('on');
      state.cine = setInterval(() => {
        state.cross.z = (state.cross.z + 1) % N;
        syncSliceUI();
        updateAll();
      }, 90);
    }
  });

  $('#sliceRange').addEventListener('input', () => {
    state.cross.z = Number($('#sliceRange').value);
    $('#sliceNum').textContent = state.cross.z + 1;
    updateAll();
  });
  $('#btnPrev').addEventListener('click', () => {
    state.cross.z = clamp(state.cross.z - 1, 0, N - 1);
    syncSliceUI(); updateAll();
  });
  $('#btnNext').addEventListener('click', () => {
    state.cross.z = clamp(state.cross.z + 1, 0, N - 1);
    syncSliceUI(); updateAll();
  });

  let rt = null;
  window.addEventListener('resize', () => {
    clearTimeout(rt);
    rt = setTimeout(() => {
      updateAll();
      render3D(true);
    }, 180);
  });

  $('#btnDemoRoi').addEventListener('click', demoRoi);
}

/* ============================================================
   AI 辅助阅片
   离线：按影像组学特征生成规则化描述
   在线：OpenAI 兼容多模态接口，发送当前轴位切片 + 特征摘要
   ============================================================ */
const AI_PRESETS = {
  gitee: { base: 'https://ai.gitee.com/v1', model: 'qwen3-vl-plus' },
  qwen: { base: 'https://dashscope.aliyuncs.com/compatible-mode/v1', model: 'qwen-vl-max' },
  zhipu: { base: 'https://open.bigmodel.cn/api/paas/v4', model: 'glm-4v-plus' },
  ark: { base: 'https://ark.cn-beijing.volces.com/api/v3', model: 'doubao-1-5-vision-pro-32k' },
  openai: { base: 'https://api.openai.com/v1', model: 'gpt-4o-mini' },
  ollama: { base: 'http://localhost:11434/v1', model: 'qwen2.5vl:7b' },
  custom: { base: '', model: '' },
};

const AI_SYSTEM = [
  '你是一名医学影像科主治医师，正在带教一名医学生阅片。请用中文、简洁、结构化地回答。',
  '输出必须包含三个小节，标题用【】包裹并单独成行：',
  '【影像所见】只描述图像中客观可观察的征象：位置、形态、边界、密度及与周围结构的关系。',
  '【分析提示】给出鉴别思路，并说明还需要哪些信息才能进一步判断。',
  '【局限性】明确指出这是 Shepp–Logan 仿真体模数据，不是真实病例，结论不可用于任何临床用途。',
  '不确定的地方直接写“不确定”，禁止编造征象或具体疾病诊断。',
].join('\n');

const ai = {
  cfg: { preset: 'qwen', base: AI_PRESETS.qwen.base, model: AI_PRESETS.qwen.model, key: '' },
  history: [],
  busy: false,
  opened: false,
};

function loadAICfg() {
  try {
    const raw = localStorage.getItem('medscope.ai');
    if (raw) Object.assign(ai.cfg, JSON.parse(raw));
  } catch (e) { /* 忽略损坏的配置 */ }
}

function saveAICfg() {
  try { localStorage.setItem('medscope.ai', JSON.stringify(ai.cfg)); } catch (e) { /* 隐私模式 */ }
}

function isLocalEndpoint(base) {
  return /^https?:\/\/(localhost|127\.0\.0\.1|\[::1\])(:|\/|$)/i.test(String(base || ''));
}

/** 本地接口（Ollama 等）允许不填 Key */
function aiOnline() {
  const c = ai.cfg;
  return Boolean(c.base && c.model && (c.key || isLocalEndpoint(c.base)));
}

/* ---------- 文本格式化 ---------- */
function escHtml(s) {
  return String(s).replace(/[&<>]/g, (m) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' }[m]));
}
function mdLite(s) {
  return escHtml(s)
    .replace(/\*\*(.+?)\*\*/g, '<b>$1</b>')
    .replace(/^【(.+?)】$/gm, '<span class="sec">$1</span>');
}

/* ---------- 统计 ---------- */
/** 层面统计：脑实质（−50~150 HU）与整体最高密度分开统计，避免颅骨污染均值 */
function sliceStats() {
  const z = state.cross.z;
  let n = 0, sum = 0, sum2 = 0, mn = Infinity, mx = -Infinity;
  let headMax = -Infinity;
  for (let y = 0; y < N; y++) {
    const base = (z * N + y) * N;
    for (let x = 0; x < N; x++) {
      const hu = huAt(base + x);
      if (hu > headMax) headMax = hu;
      if (hu < -50 || hu > 150) continue;   // 仅统计脑实质密度区间
      n++; sum += hu; sum2 += hu * hu;
      if (hu < mn) mn = hu;
      if (hu > mx) mx = hu;
    }
  }
  if (n < 20) return { n, empty: true, headMax };
  const mean = sum / n;
  const std = Math.sqrt(Math.max(0, sum2 / n - mean * mean));
  return { n, mean, std, min: mn, max: mx, headMax, empty: false };
}

function aiContext() {
  const c = state.cross;
  const st = roiStats();
  const sl = sliceStats();
  const L = [];
  L.push(`模态：CT；轴位第 ${c.z + 1}/${N} 层；窗宽 ${state.ww} HU，窗位 ${state.wl} HU；十字线 x=${c.x + 1}, y=${c.y + 1}。`);
  L.push('图像方位：按放射学惯例，图像左侧=患者右侧(R)，图像右侧=患者左侧(L)，图像上方=前(A)，下方=后(P)。');
  if (st) {
    const nx = ((state.roi.cx - (N - 1) / 2) / ((N - 1) / 2)) * 100;
    const ny = ((state.roi.cy - (N - 1) / 2) / ((N - 1) / 2)) * 100;
    L.push(`ROI 中心位于图像中心${Math.abs(nx) > 3 ? (nx > 0 ? '右' : '左') + '侧 ' + Math.abs(nx).toFixed(0) + '%' : ''}${Math.abs(ny) > 3 ? (ny < 0 ? '偏上 ' : '偏下 ') + Math.abs(ny).toFixed(0) + '%' : ''}处。`);
  }
  if (sl && !sl.empty) {
    L.push(`脑实质统计（−50~150 HU）：CT 值 ${sl.min.toFixed(0)}–${sl.max.toFixed(0)} HU，均值 ${sl.mean.toFixed(0)} HU，标准差 ${sl.std.toFixed(0)} HU；层面最高密度 ${sl.headMax.toFixed(0)} HU。`);
  } else {
    L.push('当前层面脑实质体素过少，统计从略。');
  }
  if (st) {
    L.push(`ROI：中心 (${state.roi.cx.toFixed(0)}, ${state.roi.cy.toFixed(0)})，半径 ${state.roi.r.toFixed(1)} mm，面积 ${st.n} mm²。`);
    L.push(`ROI 一阶特征：均值 ${st.mean.toFixed(1)} HU，标准差 ${st.std.toFixed(1)} HU，中位数 ${st.median.toFixed(1)} HU，P10–P90 ${st.p10.toFixed(0)}–${st.p90.toFixed(0)} HU，偏度 ${st.skew.toFixed(2)}，峰度 ${st.kurt.toFixed(2)}，熵 ${st.entropy.toFixed(2)}。`);
  } else {
    L.push('ROI：未放置。');
  }
  return L.join('\n');
}

/* ---------- 离线规则化报告 ---------- */
function offlineReport(kind) {
  const c = state.cross;
  const st = roiStats();
  const sl = sliceStats();
  if (!sl || sl.empty) {
    return '当前层面的脑实质体素过少（可能位于颅顶或颅底），统计参考价值有限。请先滚动到含脑实质的层面，或重新放置 ROI。';
  }

  if (kind === 'advice') {
    const tips = [
      '调整窗宽窗位复核边界：可试 WW 350 / WL 40（软组织窗）与 WW 2000 / WL 500（骨窗）。',
      '上下滚动相邻层面，追踪该区域的纵向范围与形态变化。',
      '若存在 ROI，可在病灶边缘与中心各取一个 ROI，比较均值与标准差。',
      '开启「体绘制」或「MIP」观察整体分布，确认是否为单发。',
    ];
    return [
      '【建议】',
      `当前层面：轴位第 ${c.z + 1} 层，窗宽 ${state.ww} / 窗位 ${state.wl} HU。`,
      ...tips.map((t, i) => `${i + 1}. ${t}`),
      '',
      '【说明】离线模式的建议由固定规则生成，未进行图像理解。',
    ].join('\n');
  }

  if (kind === 'free') {
    return [
      '离线模式下我只能按固定规则描述当前层面与 ROI 的统计特征，无法理解你提出的具体问题。',
      '',
      '可以：① 点上方快捷项获取规则化报告；② 在「设置」中填入多模态接口（Key / 接口地址 / 模型）后切换为在线识图对话。',
    ].join('\n');
  }

  const lines = [];
  lines.push('【影像信息】');
  lines.push(`CT · 轴位第 ${c.z + 1}/${N} 层 · 窗宽 ${state.ww} / 窗位 ${state.wl} HU`);
  lines.push(`脑实质体素 ${sl.n.toLocaleString('en-US')} 个（−50 ~ 150 HU）`);
  lines.push('');
  lines.push('【影像所见】');
  if (sl.headMax > 500) {
    lines.push(`可见闭合的环形高密度影，最高 ${sl.headMax.toFixed(0)} HU，符合骨性结构密度。`);
  }
  const cv = sl.mean ? sl.std / Math.abs(sl.mean) : 0;
  lines.push(`脑实质 CT 值 ${sl.min.toFixed(0)}–${sl.max.toFixed(0)} HU，均值 ${sl.mean.toFixed(0)} HU，标准差 ${sl.std.toFixed(0)} HU（变异系数 ${cv.toFixed(2)}）。`);
  if (st) {
    lines.push('');
    lines.push('【ROI 特征】');
    lines.push(`中心 (${state.roi.cx.toFixed(0)}, ${state.roi.cy.toFixed(0)})，面积 ${st.n} mm²`);
    lines.push(`均值 ${st.mean.toFixed(1)} HU，标准差 ${st.std.toFixed(1)} HU，中位数 ${st.median.toFixed(1)} HU`);
    lines.push(`P10–P90 ${st.p10.toFixed(0)}–${st.p90.toFixed(0)} HU，偏度 ${st.skew.toFixed(2)}，峰度 ${st.kurt.toFixed(2)}，熵 ${st.entropy.toFixed(2)}`);
    lines.push('');
    lines.push('【判读】');
    const delta = st.mean - sl.mean;
    const high = st.p90 - sl.mean;
    if (delta > 8) lines.push(`ROI 平均密度较周围脑实质高 ${delta.toFixed(0)} HU，存在高密度成分。`);
    else if (delta < -8) lines.push(`ROI 平均密度较周围脑实质低 ${Math.abs(delta).toFixed(0)} HU，存在低密度成分。`);
    else lines.push('ROI 平均密度与周围脑实质接近。');
    if (st.kurt < -1) lines.push(`峰度 ${st.kurt.toFixed(2)}（明显负值），密度分布呈双峰形态，提示 ROI 内包含两类不同密度的组织成分。`);
    else if (st.kurt > 1) lines.push(`峰度 ${st.kurt.toFixed(2)}（明显正值），分布尖峰，说明局部密度高度集中。`);
    if (Math.abs(st.skew) > 0.8) lines.push(`偏度 ${st.skew.toFixed(2)}，分布${st.skew > 0 ? '右偏（高密度侧拖尾）' : '左偏（低密度侧拖尾）'}。`);
    if (high > 30) lines.push(`P90 高出脑实质均值 ${high.toFixed(0)} HU，高密度区占 ROI 的比例不可忽略。`);
  } else {
    lines.push('');
    lines.push('【ROI 特征】');
    lines.push('尚未放置 ROI。切换到「圆形 ROI」工具，在轴位视图上拖动画圆即可统计。');
  }
  lines.push('');
  lines.push('【说明】');
  lines.push('以上内容由固定规则依据当前窗宽窗位与 ROI 统计生成（离线模式），未经图像理解，不构成诊断。');
  return lines.join('\n');
}

/* ---------- 消息渲染 ---------- */
function aiScroll() {
  const body = $('#aiBody');
  body.scrollTop = body.scrollHeight;
}

function aiPush(role, html, cls = '') {
  const wrap = document.createElement('div');
  wrap.className = 'msg ' + role;
  const bubble = document.createElement('div');
  bubble.className = 'bubble' + (cls ? ' ' + cls : '');
  if (role === 'user') bubble.textContent = html;
  else bubble.innerHTML = html;
  wrap.appendChild(bubble);
  $('#aiBody').appendChild(wrap);
  aiScroll();
  return bubble;
}

function aiWelcome() {
  const online = aiOnline();
  aiPush('bot', mdLite([
    '【使用方式】',
    online
      ? `已连接多模态接口：${escHtml(ai.cfg.model)}。点下方快捷项或直接追问，我会读取当前轴位切片并作答。`
      : '当前为离线模式：我会基于当前窗宽窗位、ROI 的影像组学特征给出规则化描述；如需真正的识图问答，请在「设置」中填入多模态接口。',
    '',
    '【可问什么】',
    '· 描述当前层面可见的征象',
    '· 选中病灶区后询问其密度特征',
    '· 下一步应补充哪些观察或检查',
    '',
    '【提示】',
    '先用「圆形 ROI」在病灶上画一个圈，提问会更有针对性。',
  ].join('\n')));
}

/* ---------- 发送 ---------- */
async function aiAsk(text, kind) {
  if (ai.busy) return;
  const question = text.trim();
  if (!question) return;
  aiPush('user', question);
  ai.history.push({ role: 'user', content: question });
  $('#aiText').value = '';

  if (!aiOnline()) {
    const reply = offlineReport(kind || 'free');
    ai.history.push({ role: 'assistant', content: reply });
    setTimeout(() => aiPush('bot', mdLite(reply)), 120);
    return;
  }

  ai.busy = true;
  $('#aiSend').disabled = true;
  const bubble = aiPush('bot', '<span style="color:#5d6e7e">正在分析…</span>', 'typing');

  try {
    const img = captureSlice('axial', 512);
    const msgs = [{ role: 'system', content: AI_SYSTEM }];
    for (const h of ai.history.slice(-8)) msgs.push({ role: h.role, content: h.content });
    msgs.push({
      role: 'user',
      content: [
        { type: 'text', text: `${aiContext()}\n\n问题：${question}` },
        { type: 'image_url', image_url: { url: img } },
      ],
    });

    const headers = { 'Content-Type': 'application/json' };
    if (ai.cfg.key) headers.Authorization = 'Bearer ' + ai.cfg.key;

    const res = await fetch(ai.cfg.base.replace(/\/+$/, '') + '/chat/completions', {
      method: 'POST',
      headers,
      body: JSON.stringify({ model: ai.cfg.model, stream: true, max_tokens: 900, temperature: 0.3, messages: msgs }),
    });

    if (!res.ok) {
      const detail = await res.text().catch(() => '');
      throw new Error(`HTTP ${res.status} ${detail.slice(0, 160)}`);
    }

    bubble.classList.remove('typing');
    bubble.innerHTML = '';
    let acc = '';
    const reader = res.body.getReader();
    const dec = new TextDecoder();
    let buf = '';
    for (;;) {
      const { value, done } = await reader.read();
      if (done) break;
      buf += dec.decode(value, { stream: true });
      let idx;
      while ((idx = buf.indexOf('\n')) >= 0) {
        const line = buf.slice(0, idx).trim();
        buf = buf.slice(idx + 1);
        if (!line.startsWith('data:')) continue;
        const payload = line.slice(5).trim();
        if (payload === '[DONE]') continue;
        try {
          const j = JSON.parse(payload);
          const delta = j.choices?.[0]?.delta?.content;
          if (delta) {
            acc += delta;
            bubble.textContent = acc;
            aiScroll();
          }
        } catch (e) { /* 忽略心跳行 */ }
      }
    }
    bubble.innerHTML = mdLite(acc || '（接口未返回内容）');
    ai.history.push({ role: 'assistant', content: acc });
    aiScroll();
  } catch (err) {
    bubble.classList.remove('typing');
    bubble.classList.add('err');
    const msg = String(err && err.message ? err.message : err);
    let hint = '';
    if (/Failed to fetch|NetworkError|Load failed/i.test(msg)) {
      hint = '可能是跨域被拒绝：直接用 file:// 打开时浏览器会屏蔽该请求。可用 python -m http.server 起本地服务后访问 http://localhost:8000/02-Imaging/index.html，或换用允许浏览器调用的接口。';
    } else if (/401|403/.test(msg)) {
      hint = 'API Key 无效或没有该模型的权限。';
    } else if (/404/.test(msg)) {
      hint = '接口地址或模型名称不正确。';
    } else if (/429/.test(msg)) {
      hint = '请求过于频繁或余额不足。';
    }
    bubble.innerHTML = mdLite(`【请求失败】\n${msg}${hint ? '\n\n【排查建议】\n' + hint : ''}`);
    aiScroll();
  } finally {
    ai.busy = false;
    $('#aiSend').disabled = false;
  }
}

/* ---------- 切片截图 ---------- */
function captureSlice(view, size) {
  fillSlice(view);
  const c = document.createElement('canvas');
  c.width = size;
  c.height = size;
  const x = c.getContext('2d');
  x.imageSmoothingEnabled = true;
  x.imageSmoothingQuality = 'high';
  x.drawImage(sliceCanvas, 0, 0, size, size);
  return c.toDataURL('image/png');
}

/* ---------- 设置 ---------- */
function syncAIConfigUI() {
  $('#cfgPreset').value = ai.cfg.preset || 'qwen';
  $('#cfgBase').value = ai.cfg.base || '';
  $('#cfgModel').value = ai.cfg.model || '';
  $('#cfgKey').value = ai.cfg.key || '';
}

function syncAIStatus() {
  const on = aiOnline();
  const local = isLocalEndpoint(ai.cfg.base);
  $('#aiState').classList.toggle('online', on);
  $('#aiDot').style.background = on ? 'var(--teal)' : 'var(--cyan)';
  $('#aiMode').textContent = on
    ? `${local ? '本地' : '在线'} · ${ai.cfg.model}`
    : '离线模式 · 规则化描述';
}

/* ---------- 连接测试 ---------- */
async function aiTestConnection() {
  const out = $('#cfgTestResult');
  const base = $('#cfgBase').value.trim().replace(/\/+$/, '');
  const model = $('#cfgModel').value.trim();
  const key = $('#cfgKey').value.trim();
  if (!base || !model) {
    out.className = 'test-result bad';
    out.textContent = '请先填写接口地址与模型名称';
    return;
  }
  out.className = 'test-result';
  out.textContent = '测试中…';
  const btn = $('#cfgTest');
  btn.disabled = true;
  try {
    const headers = { 'Content-Type': 'application/json' };
    if (key) headers.Authorization = 'Bearer ' + key;
    const res = await fetch(base + '/chat/completions', {
      method: 'POST',
      headers,
      body: JSON.stringify({
        model,
        max_tokens: 16,
        stream: false,
        messages: [{ role: 'user', content: '回复两个字：可用' }],
      }),
    });
    if (!res.ok) {
      const t = await res.text().catch(() => '');
      throw new Error(`HTTP ${res.status} ${t.slice(0, 120)}`);
    }
    const j = await res.json();
    const txt = (j.choices?.[0]?.message?.content || '').trim().slice(0, 40);
    out.className = 'test-result ok';
    out.textContent = `连接成功${txt ? '：' + txt : ''}`;
  } catch (err) {
    const msg = String(err && err.message ? err.message : err);
    out.className = 'test-result bad';
    out.textContent = /Failed to fetch|NetworkError|Load failed/i.test(msg)
      ? '无法连接（可能是跨域、地址错误或服务未启动）'
      : msg.slice(0, 90);
  } finally {
    btn.disabled = false;
  }
}

function bindAI() {
  loadAICfg();
  syncAIStatus();

  $('#aiOpen').addEventListener('click', () => {
    const d = $('#aiDrawer');
    const open = !d.classList.contains('open');
    d.classList.toggle('open', open);
    d.setAttribute('aria-hidden', String(!open));
    $('#aiOpen').classList.toggle('on', open);
    if (open && !ai.opened) {
      ai.opened = true;
      aiWelcome();
    }
  });
  $('#aiClose').addEventListener('click', () => {
    $('#aiDrawer').classList.remove('open');
    $('#aiDrawer').setAttribute('aria-hidden', 'true');
    $('#aiOpen').classList.remove('on');
  });

  $('#aiForm').addEventListener('submit', (e) => {
    e.preventDefault();
    aiAsk($('#aiText').value, 'free');
  });
  $('#aiText').addEventListener('keydown', (e) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      aiAsk($('#aiText').value, 'free');
    }
  });

  $$('#aiQuick button').forEach((b) => {
    b.addEventListener('click', () => {
      const kind = b.dataset.q;
      const map = {
        describe: '请描述当前层面的影像所见。',
        roi: '请分析当前 ROI 的密度特征，并说明可能提示什么。',
        advice: '针对当前所见，下一步应该补充哪些观察或检查？',
      };
      aiAsk(map[kind], kind);
    });
  });

  // 设置弹窗
  const modal = $('#aiModal');
  const openModal = () => { syncAIConfigUI(); modal.classList.add('open'); modal.setAttribute('aria-hidden', 'false'); };
  const closeModal = () => { modal.classList.remove('open'); modal.setAttribute('aria-hidden', 'true'); };

  $('#aiSettings').addEventListener('click', openModal);
  $('#cfgCancel').addEventListener('click', closeModal);
  $('#cfgTest').addEventListener('click', aiTestConnection);
  modal.addEventListener('click', (e) => { if (e.target === modal) closeModal(); });

  $('#cfgPreset').addEventListener('change', () => {
    const p = AI_PRESETS[$('#cfgPreset').value];
    if (p && $('#cfgPreset').value !== 'custom') {
      $('#cfgBase').value = p.base;
      $('#cfgModel').value = p.model;
    }
  });

  $('#cfgSave').addEventListener('click', () => {
    ai.cfg.preset = $('#cfgPreset').value;
    ai.cfg.base = $('#cfgBase').value.trim();
    ai.cfg.model = $('#cfgModel').value.trim();
    ai.cfg.key = $('#cfgKey').value.trim();
    saveAICfg();
    syncAIStatus();
    closeModal();
    aiPush('bot', mdLite(aiOnline()
      ? `【已切换为在线模式】接口：${escHtml(ai.cfg.base)}，模型：${escHtml(ai.cfg.model)}。`
      : '【已切换为离线模式】未填写 API Key，将继续使用规则化描述。'));
  });

  $('#cfgClear').addEventListener('click', () => {
    ai.cfg = { preset: 'qwen', base: AI_PRESETS.qwen.base, model: AI_PRESETS.qwen.model, key: '' };
    saveAICfg();
    syncAIConfigUI();
    syncAIStatus();
    closeModal();
    aiPush('bot', '【已清除配置】当前为离线模式。');
  });

  document.addEventListener('keydown', (e) => {
    if (e.key === 'Escape') {
      if (modal.classList.contains('open')) closeModal();
      else if ($('#aiDrawer').classList.contains('open')) $('#aiClose').click();
    }
  });
}

/* ---------------- 启动 ---------------- */
(function boot() {
  const overlay = document.createElement('div');
  overlay.className = 'boot';
  overlay.innerHTML = '<div class="boot-inner"><div class="boot-ring"></div><p>正在重建体数据…</p></div>';
  document.body.appendChild(overlay);

  setTimeout(() => {
    const t0 = performance.now();
    buildVolume();
    buildLUT();
    const dt = performance.now() - t0;

    bindControls();
    bindViewports();
    bindAI();
    syncSliceUI();
    buildFilm();
    updateAll();
    setTimeout(() => render3D(true), 60);

    overlay.classList.add('done');
    setTimeout(() => overlay.remove(), 400);
    console.log(`[MedScope] volume ${N}³ built in ${dt.toFixed(0)} ms`);

    if (location.hash === '#roi') demoRoi();
  }, 40);
})();

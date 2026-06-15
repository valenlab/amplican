/*
  amplican_animation.js
  Self-contained SVG animation that walks through the amplican pipeline.
  Pairs with inst/amplican_animation.html. No external dependencies.
  Data is taken from the package's own example dataset (inst/extdata),
  primarily experiment ID_1, so every counter matches config_summary.csv.
*/
(function () {
"use strict";

/* ------------------------------------------------------------------ */
/* Engine                                                              */
/* ------------------------------------------------------------------ */
const SVGNS = "http://www.w3.org/2000/svg";
const canvas = document.getElementById("canvas");
const canvasWrap = document.getElementById("canvasWrap");

const lerp = (a, b, t) => a + (b - a) * t;
const easeOut = t => 1 - Math.pow(1 - t, 3);
const easeInOut = t => (t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2);
const easeBack = t => {
  const c1 = 1.70158, c3 = c1 + 1;
  return 1 + c3 * Math.pow(t - 1, 3) + c1 * Math.pow(t - 1, 2);
};

let token = { cancelled: false };

function E(tag, attrs, parent) {
  const el = document.createElementNS(SVGNS, tag);
  el.state = { x: 0, y: 0, scale: 1, opacity: 1 };
  if (attrs) for (const k in attrs) {
    if (k === "text") el.textContent = attrs[k];
    else el.setAttribute(k, attrs[k]);
  }
  if (parent) parent.appendChild(el);
  return el;
}
function place(el, x, y, scale) {
  const s = el.state;
  if (x !== undefined) s.x = x;
  if (y !== undefined) s.y = y;
  if (scale !== undefined) s.scale = scale;
  el.setAttribute("transform", "translate(" + s.x + " " + s.y + ") scale(" + s.scale + ")");
  return el;
}
function setOpacity(el, o) { el.style.opacity = o; el.state.opacity = o; return el; }

function applyAnimProp(el, k) {
  const s = el.state;
  if (k === "x" || k === "y" || k === "scale")
    el.setAttribute("transform", "translate(" + s.x + " " + s.y + ") scale(" + s.scale + ")");
  else if (k === "opacity") el.style.opacity = s.opacity;
  else el.setAttribute(k, s[k]);
}
function anim(el, target, duration, easing) {
  if (duration == null) duration = 600;
  if (easing == null) easing = easeOut;
  const from = {};
  for (const k in target) from[k] = (el.state[k] != null ? el.state[k] : 0);
  const myToken = token;
  return new Promise(function (resolve) {
    const t0 = performance.now();
    function frame(now) {
      if (myToken.cancelled) {
        for (const k in target) { el.state[k] = target[k]; applyAnimProp(el, k); }
        resolve(); return;
      }
      const t = Math.min(1, (now - t0) / duration);
      const e = easing(t);
      for (const k in target) {
        el.state[k] = from[k] + (target[k] - from[k]) * e;
        applyAnimProp(el, k);
      }
      if (t < 1) requestAnimationFrame(frame); else resolve();
    }
    requestAnimationFrame(frame);
  });
}
function sleep(ms) {
  const myToken = token;
  return new Promise(function (resolve) {
    const t0 = performance.now();
    function check(now) {
      if (myToken.cancelled) { resolve(); return; }
      if (now - t0 >= ms) resolve(); else requestAnimationFrame(check);
    }
    requestAnimationFrame(check);
  });
}
function fadeIn(el, dur) { setOpacity(el, 0); return anim(el, { opacity: 1 }, dur || 400); }
function fadeOut(el, dur) { return anim(el, { opacity: 0 }, dur || 400); }
function popIn(el, x, y, dur) {
  setOpacity(el, 0); place(el, x, y, 0.6);
  return anim(el, { opacity: 1, scale: 1 }, dur || 500, easeBack);
}

/* ------------------------------------------------------------------ */
/* Palette + layout                                                    */
/* ------------------------------------------------------------------ */
const PAL = {
  ink: "#0f172a", muted: "#64748b", line: "#e2e8f0", soft: "#f1f5f9",
  primer: "#10b981", cut: "#f59e0b", cutwin: "#fde68a",
  del: "#ef4444", ins: "#8b5cf6", hdr: "#0ea5e9", pink: "#ec4899",
  A: "#ef4444", C: "#3b82f6", G: "#eab308", T: "#22c55e", N: "#94a3b8"
};
const nucColor = b => PAL[b] || PAL.N;

const AX0 = 70, AX1 = 930;       // amplicon horizontal extent (viewBox coords)
const AMPLEN = 199;
const bp2x = bp => AX0 + (bp - 1) / (AMPLEN - 1) * (AX1 - AX0);

/* ------------------------------------------------------------------ */
/* Data (from inst/extdata)                                            */
/* ------------------------------------------------------------------ */
const AMPLICON = {
  id: "ID_1", length: 199,
  fwdPrimer: "AAGCTGACGGCTAAATGA", fwd: [1, 18],
  revPrimerRC: "GTGTGTTTGCGCTTGTGTAATT", rev: [178, 199],
  guide: "AGGTGGTCAGGGAACTGG", cut: [66, 83],
  buffer: 5, window: [61, 88]
};

// Six representative "fates" used in the flow scenes.
const CAST = {
  R1: { id: "R1", label: "low quality", count: 3, fill: "#fecaca", stroke: "#ef4444", fate: "bad average quality" },
  R2: { id: "R2", label: "unassigned", count: 4, fill: "#fed7aa", stroke: "#f97316", fate: "no primer match" },
  R3: { id: "R3", label: "primer dimer", count: 1, fill: "#ede9fe", stroke: "#8b5cf6", fate: "primer dimer",
        events: [{ type: "deletion", start: 24, end: 173, width: 150 }] },
  R4: { id: "R4", label: "wild-type", count: 1, fill: "#dcfce7", stroke: "#22c55e", fate: "no edits" },
  R5: { id: "R5", label: "edited", count: 3, fill: "#fce7f3", stroke: "#ec4899", fate: "deletion at cut",
        events: [{ type: "deletion", start: 42, end: 107, width: 66 },
                 { type: "insertion", start: 108, end: 127, width: 20, alt: "A" },
                 { type: "mismatch", start: 31, end: 31, ref: "A", alt: "G" },
                 { type: "mismatch", start: 163, end: 163, ref: "T", alt: "A" }] },
  R6: { id: "R6", label: "HDR", count: 2, fill: "#e0f2fe", stroke: "#0ea5e9", fate: "HDR",
        events: [{ type: "mismatch", start: 121, end: 121, ref: "G", alt: "A" }], hdr: true }
};

// Full ID_1 read set (forward strand) used in alignment / summarization scenes.
const ID1_READS = [
  { id: "read 1", count: 3, score: 597, events: CAST.R5.events },
  { id: "read 2", count: 2, score: 557,
    events: [{ type: "deletion", start: 38, end: 111, width: 74 },
             { type: "insertion", start: 112, end: 131, width: 20, alt: "A" },
             { type: "mismatch", start: 31, end: 31, ref: "A", alt: "G" },
             { type: "mismatch", start: 163, end: 163, ref: "T", alt: "A" }] },
  { id: "read 3", count: 1, score: 193, events: CAST.R3.events },        // primer dimer
  { id: "read 4", count: 4, countOrig: 1, score: 532,
    events: [{ type: "deletion", start: 34, end: 117, width: 84 },
             { type: "mismatch", start: 31, end: 31, ref: "A", alt: "G" },
             { type: "mismatch", start: 171, end: 171, ref: "G", alt: "A" }] }
];
// Note: read counts above are the real unique-read weights from raw_events.csv
// (3 + 2 + 1 + 1 = 7 = ID_1 Reads).

const ID1_SUMMARY = {
  Reads: 7, PRIMER_DIMER: 1, Low_Score: 0, Reads_Filtered: 6,
  Reads_Del: 6, Reads_In: 0, Reads_Edited: 6, Reads_Frameshifted: 2, HDR: 0
};

/* ------------------------------------------------------------------ */
/* Small drawing helpers                                               */
/* ------------------------------------------------------------------ */
function label(x, y, text, opts) {
  opts = opts || {};
  const t = E("text", {
    x: x, y: y, "text-anchor": opts.anchor || "start",
    "font-size": opts.size || 12, "font-weight": opts.weight || 400,
    fill: opts.fill || PAL.muted
  }, opts.parent || canvas);
  t.textContent = text;
  if (opts.opacity != null) setOpacity(t, opts.opacity);
  return t;
}

function arrow(x1, y1, x2, y2, color, marker) {
  return E("line", {
    x1: x1, y1: y1, x2: x2, y2: y2,
    stroke: color || "#94a3b8", "stroke-width": 1.6,
    "marker-end": "url(#" + (marker || "arr") + ")"
  }, canvas);
}

function Counter(x, y, labelTxt, value, color) {
  const grp = E("g", {}, canvas); place(grp, x, y);
  E("rect", { x: -52, y: -24, width: 104, height: 48, rx: 9, fill: "#fff",
    stroke: color || PAL.line, "stroke-width": 1.5 }, grp);
  const valEl = E("text", { x: 0, y: 4, "text-anchor": "middle", "font-size": 22,
    "font-weight": 700, fill: color || PAL.ink }, grp);
  valEl.textContent = value;
  const labEl = E("text", { x: 0, y: 22, "text-anchor": "middle", "font-size": 9.5,
    fill: PAL.muted }, grp);
  labEl.textContent = labelTxt;
  const o = { grp: grp, value: value, set };
  function set(v, dur) {
    const from = o.value;
    const myToken = token;
    return new Promise(function (resolve) {
      const t0 = performance.now();
      function f(now) {
        if (myToken.cancelled) { o.value = v; valEl.textContent = v; resolve(); return; }
        const t = Math.min(1, (now - t0) / (dur || 500));
        o.value = Math.round(from + (v - from) * t);
        valEl.textContent = o.value;
        if (t < 1) requestAnimationFrame(f); else resolve();
      }
      requestAnimationFrame(f);
    });
  }
  return o;
}

function Bin(x, y, w, h, labelTxt, stroke) {
  const grp = E("g", {}, canvas); place(grp, x, y);
  E("rect", { x: 0, y: 0, width: w, height: h, rx: 12, fill: "#fff",
    stroke: stroke || "#cbd5e1", "stroke-width": 1.4, "stroke-dasharray": "5 4" }, grp);
  const t = E("text", { x: w / 2, y: h + 16, "text-anchor": "middle", "font-size": 11,
    fill: PAL.muted }, grp);
  t.textContent = labelTxt;
  return grp;
}

// A compact read token (not bp-aligned) for flow / conveyor scenes.
function readToken(read, opts) {
  opts = opts || {};
  const w = opts.w || 78;
  const grp = E("g", {}, canvas);
  E("rect", { x: -w / 2, y: -9, width: w, height: 18, rx: 9,
    fill: read.fill || PAL.soft, stroke: read.stroke || "#94a3b8", "stroke-width": 1.2 }, grp);
  const t = E("text", { x: 0, y: 3, "text-anchor": "middle", "font-size": 10,
    "font-weight": 600, fill: PAL.ink }, grp);
  t.textContent = read.id + (read.count > 1 ? "  ×" + read.count : "");
  place(grp, opts.x || 0, opts.y || 0);
  if (opts.opacity != null) setOpacity(grp, opts.opacity);
  return grp;
}

// Amplicon track. Returns handles to the regions so scenes can highlight them.
function drawAmplicon(y, opts) {
  opts = opts || {};
  const g = E("g", {}, canvas); place(g, 0, 0);
  const cy = y;
  // backbone
  E("rect", { x: AX0, y: cy - 7, width: AX1 - AX0, height: 14, rx: 7,
    fill: PAL.soft, stroke: PAL.line }, g);
  // cut window band
  if (opts.window !== false) {
    E("rect", { x: bp2x(AMPLICON.window[0]), y: cy - 11,
      width: bp2x(AMPLICON.window[1]) - bp2x(AMPLICON.window[0]), height: 22,
      rx: 4, fill: PAL.cutwin, opacity: 0.7 }, g);
  }
  // forward primer
  E("rect", { x: bp2x(AMPLICON.fwd[0]), y: cy - 7,
    width: bp2x(AMPLICON.fwd[1]) - bp2x(AMPLICON.fwd[0]), height: 14, rx: 3, fill: PAL.primer }, g);
  // reverse primer RC
  E("rect", { x: bp2x(AMPLICON.rev[0]), y: cy - 7,
    width: bp2x(AMPLICON.rev[1]) - bp2x(AMPLICON.rev[0]) + 2, height: 14, rx: 3, fill: PAL.primer }, g);
  // cut site / guide (UPPERCASE)
  const cutRect = E("rect", { x: bp2x(AMPLICON.cut[0]), y: cy - 7,
    width: bp2x(AMPLICON.cut[1]) - bp2x(AMPLICON.cut[0]), height: 14, rx: 3, fill: PAL.cut }, g);
  // position ticks (boundaries of the primer / guide regions)
  [AMPLICON.fwd[1], AMPLICON.cut[0], AMPLICON.cut[1], AMPLICON.rev[0]].forEach(function (bp) {
    E("line", { x1: bp2x(bp), y1: cy + 8, x2: bp2x(bp), y2: cy + 14, stroke: PAL.muted }, g);
  });
  return { g: g, cut: cutRect, y: cy };
}

function annotateAmplicon(amp, showLabels) {
  if (!showLabels) return;
  label(bp2x((AMPLICON.fwd[0] + AMPLICON.fwd[1]) / 2), amp.y - 16, "Forward primer",
    { anchor: "middle", size: 10.5, fill: "#047857", weight: 600 });
  label(bp2x((AMPLICON.rev[0] + AMPLICON.rev[1]) / 2), amp.y - 16, "Reverse primer (RC)",
    { anchor: "middle", size: 10.5, fill: "#047857", weight: 600 });
  label(bp2x((AMPLICON.cut[0] + AMPLICON.cut[1]) / 2), amp.y + 30, "guide / cut site (UPPERCASE)",
    { anchor: "middle", size: 10.5, fill: "#b45309", weight: 600 });
  label(bp2x(1), amp.y + 30, "1", { anchor: "middle", size: 9.5 });
  label(bp2x(AMPLICON.length), amp.y + 30, String(AMPLICON.length), { anchor: "middle", size: 9.5 });
}

// Subtract deletions from [1,len] -> surviving segments.
function nonDeleted(len, dels) {
  let ints = [[1, len]];
  (dels || []).forEach(function (d) {
    const s = d.start, e = d.end;
    const next = [];
    ints.forEach(function (it) {
      if (it[1] < s || it[0] > e) next.push(it);
      else {
        if (it[0] < s) next.push([it[0], s - 1]);
        if (it[1] > e) next.push([e + 1, it[1]]);
      }
    });
    ints = next;
  });
  return ints;
}

// Detailed, bp-aligned read under the amplicon.
function detailedRead(read, y, opts) {
  opts = opts || {};
  const g = E("g", {}, canvas); place(g, 0, 0);
  const dels = (read.events || []).filter(e => e.type === "deletion");
  const segs = nonDeleted(AMPLEN, dels);
  // segments
  segs.forEach(function (it) {
    const x1 = bp2x(it[0]), x2 = bp2x(it[1]);
    E("rect", { x: x1, y: y - 5, width: Math.max(1.5, x2 - x1), height: 10,
      rx: 3, fill: opts.fill || "#cbd5e1", stroke: "#94a3b8", "stroke-width": 0.8 }, g);
  });
  // events
  (read.events || []).forEach(function (ev) {
    if (ev.type === "mismatch") {
      E("circle", { cx: bp2x(ev.start), cy: y, r: 4.5, fill: nucColor(ev.alt),
        stroke: "#fff", "stroke-width": 1 }, g);
    } else if (ev.type === "insertion") {
      const px = bp2x((ev.start + ev.end) / 2);
      E("polygon", {
        points: px + "," + (y + 5) + " " + (px - 5) + "," + (y + 13) + " " + (px + 5) + "," + (y + 13),
        fill: PAL.ins
      }, g);
    } else if (ev.type === "deletion") {
      const x1 = bp2x(ev.start), x2 = bp2x(ev.end);
      E("rect", { x: x1, y: y - 13, width: Math.max(2, x2 - x1), height: 6, rx: 2,
        fill: PAL.del, opacity: 0.85 }, g);
      E("line", { x1: x1, y1: y - 10, x2: x2, y2: y - 10, stroke: PAL.del, "stroke-width": 1,
        "stroke-dasharray": "2 2" }, g);
    }
  });
  // id label
  const t = E("text", { x: AX0 - 14, y: y + 4, "text-anchor": "end", "font-size": 11,
    fill: PAL.ink, "font-weight": 600 }, g);
  t.textContent = read.id + (read.count > 1 ? "  ×" + read.count : "");
  return g;
}

// event "chip" used in the extraction scene
function chip(x, y, text, color) {
  const g = E("g", {}, canvas);
  const pad = 6;
  // measure approx
  const w = text.length * 6.2 + pad * 2;
  E("rect", { x: 0, y: 0, width: w, height: 22, rx: 6, fill: "#fff",
    stroke: color || PAL.line, "stroke-width": 1.2 }, g);
  const t = E("text", { x: w / 2, y: 14, "text-anchor": "middle", "font-size": 10.5,
    fill: PAL.ink, "font-weight": 600 }, g);
  t.textContent = text;
  E("circle", { cx: 8, cy: 11, r: 3.5, fill: color || PAL.muted }, g);
  place(g, x - w / 2, y);
  return g;
}

/* ------------------------------------------------------------------ */
/* Scene controller                                                    */
/* ------------------------------------------------------------------ */
const figurePane = document.getElementById("figurePane");
function clearFigures() { figurePane.innerHTML = ""; }
function showFigure(src, caption, opts) {
  opts = opts || {};
  const wrap = document.createElement("div");
  wrap.className = "fig";
  wrap.style.opacity = "0";
  const img = document.createElement("img");
  img.src = src;
  img.alt = caption || "";
  wrap.appendChild(img);
  if (caption) {
    const c = document.createElement("div");
    c.className = "fcap";
    c.textContent = caption;
    wrap.appendChild(c);
  }
  figurePane.appendChild(wrap);
  // fade in
  const t0 = performance.now();
  const myToken = token;
  function frame(now) {
    if (myToken.cancelled) { wrap.style.opacity = 1; return; }
    const t = Math.min(1, (now - t0) / (opts.dur || 500));
    wrap.style.opacity = t;
    if (t < 1) requestAnimationFrame(frame);
  }
  requestAnimationFrame(frame);
  return wrap;
}

function clearStage() {
  canvas.setAttribute("viewBox", "0 0 1000 560");
  canvas.style.display = "";
  if (canvasWrap) canvasWrap.style.display = "";
  for (let i = canvas.children.length - 1; i >= 0; i--) {
    const child = canvas.children[i];
    if (child.tagName !== "defs") child.remove();
  }
  clearFigures();
}

const SCENES = [];
function defineScene(meta, run) { SCENES.push(Object.assign({ run: run }, meta)); }

// Token-aware, cancellable scheduler. Signature matches setTimeout (fn, ms)
// so staggered scene work can be cancelled on scene change and never leak
// into the next scene's canvas.
const pendingTimeouts = [];
function later(fn, ms) {
  const myToken = token;
  const id = setTimeout(function () {
    const k = pendingTimeouts.indexOf(id);
    if (k >= 0) pendingTimeouts.splice(k, 1);
    if (!myToken.cancelled) fn();
  }, ms);
  pendingTimeouts.push(id);
  return id;
}

let idx = 0;
let currentRun = null;     // promise of the scene currently playing
async function show(i) {
  i = Math.max(0, Math.min(SCENES.length - 1, i));
  idx = i;
  token.cancelled = true;          // cancel the running scene's awaits/timers
  token = { cancelled: false };
  const myToken = token;
  while (pendingTimeouts.length) clearTimeout(pendingTimeouts.pop());
  // Let the previous scene's run() unwind first: once cancelled, its awaits
  // resolve instantly so it races to completion in a few ms. Awaiting it here
  // means any straggler appends happen BEFORE we wipe the canvas — so they
  // can't leak into the next scene.
  if (currentRun) { try { await currentRun; } catch (e) { /* cancelled */ } }
  if (myToken.cancelled) return;   // another navigation overtook us
  clearStage();
  updateNav();
  await sleep(40);
  if (myToken.cancelled) return;
  currentRun = SCENES[i].run(myToken);
  try { await currentRun; }
  catch (e) { console.error("scene " + i + " failed", e); }
  currentRun = null;
}

function updateNav() {
  const s = SCENES[idx];
  document.getElementById("sceneTitle").textContent = (idx + 1) + ". " + s.title;
  document.getElementById("stepNo").textContent = "Step " + (idx + 1) + " / " + SCENES.length;
  document.getElementById("capTag").textContent = s.tag;
  document.getElementById("capTitle").textContent = s.title;
  document.getElementById("capBody").innerHTML = s.body;
  document.getElementById("progress").style.setProperty("--p", ((idx) / (SCENES.length - 1)) * 100 + "%");
  document.getElementById("prev").disabled = idx === 0;
  document.getElementById("next").disabled = idx === SCENES.length - 1;
  [].forEach.call(document.querySelectorAll("#chapters button"), function (b, i2) {
    b.classList.toggle("active", i2 === idx);
  });
}

function buildChapters() {
  const nav = document.getElementById("chapters");
  SCENES.forEach(function (s, i) {
    const b = document.createElement("button");
    b.textContent = (i + 1) + ". " + s.title;
    b.addEventListener("click", function () { show(i); });
    nav.appendChild(b);
  });
}

document.getElementById("next").addEventListener("click", () => show(idx + 1));
document.getElementById("prev").addEventListener("click", () => show(idx - 1));
document.getElementById("replay").addEventListener("click", () => show(idx));
window.addEventListener("keydown", function (e) {
  if (e.key === "ArrowRight" || e.key === "PageDown") { e.preventDefault(); show(idx + 1); }
  else if (e.key === "ArrowLeft" || e.key === "PageUp") { e.preventDefault(); show(idx - 1); }
  else if (e.key === " ") { e.preventDefault(); show(idx + 1); }
  else if (e.key === "Home") { e.preventDefault(); show(0); }
  else if (e.key === "End") { e.preventDefault(); show(SCENES.length - 1); }
  else if (e.key === "r" || e.key === "R") { e.preventDefault(); show(idx); }
});

/* ================================================================== */
/* SCENES                                                              */
/* ================================================================== */

/* 1. Read & amplicon --------------------------------------------- */
defineScene({
  tag: "setup",
  title: "The read & its amplicon",
  body: "The amplicon carries the <strong>forward primer</strong> (left), the <strong>reverse primer</strong> reverse-complemented (right), and the <strong>UPPERCASE</strong> guide / cut site in the middle. A sequencing <strong>read</strong> of the same molecule carries the same primers at its ends — that is what lets amplican anchor it later."
}, async function () {
  label(500, 40, "amplicon (ID_1)", { anchor: "middle", size: 13, weight: 700 });

  const amp = drawAmplicon(150);
  setOpacity(amp.g, 0);
  await fadeIn(amp.g, 700);
  // minimal region labels (no bp numbers)
  label(bp2x((AMPLICON.fwd[0] + AMPLICON.fwd[1]) / 2), 150 - 18, "forward primer",
    { anchor: "middle", size: 10.5, fill: "#047857", weight: 600 });
  label(bp2x((AMPLICON.rev[0] + AMPLICON.rev[1]) / 2), 150 - 18, "reverse primer (RC)",
    { anchor: "middle", size: 10.5, fill: "#047857", weight: 600 });
  label(bp2x((AMPLICON.cut[0] + AMPLICON.cut[1]) / 2), 150 + 30, "guide / cut site (UPPERCASE)",
    { anchor: "middle", size: 10.5, fill: "#b45309", weight: 600 });

  await sleep(300);
  // the read, below, with its own primers at the ends
  const ry = 290;
  // backbone
  E("rect", { x: AX0, y: ry - 6, width: AX1 - AX0, height: 12, rx: 6, fill: PAL.soft, stroke: PAL.line }, canvas);
  // read's own primers (green) at the ends — same idea, on the read
  E("rect", { x: bp2x(AMPLICON.fwd[0]), y: ry - 6, width: bp2x(AMPLICON.fwd[1]) - bp2x(AMPLICON.fwd[0]), height: 12, rx: 3, fill: PAL.primer }, canvas);
  E("rect", { x: bp2x(AMPLICON.rev[0]), y: ry - 6, width: bp2x(AMPLICON.rev[1]) - bp2x(AMPLICON.rev[0]) + 2, height: 12, rx: 3, fill: PAL.primer }, canvas);
  label(500, ry - 22, "sequencing read", { anchor: "middle", size: 12, weight: 700 });

  // connectors hinting primers match between read and amplicon
  [[AMPLICON.fwd, "left"], [AMPLICON.rev, "right"]].forEach(p => {
    const x = bp2x((p[0][0] + p[0][1]) / 2);
    E("line", { x1: x, y1: 157, x2: x, y2: ry - 7, stroke: "#cbd5e1", "stroke-width": 1, "stroke-dasharray": "3 3" }, canvas);
  });

  // the cut-site region is where editing happens — soft highlight on both
  E("rect", { x: bp2x(AMPLICON.cut[0]), y: ry - 9, width: bp2x(AMPLICON.cut[1]) - bp2x(AMPLICON.cut[0]), height: 18, rx: 3, fill: PAL.cut, opacity: 0.35 }, canvas);
});

/* 2. Quality scan ----------------------------------------------- */
defineScene({
  tag: "quality",
  title: "Quality scan",
  body: "A cursor sweeps across the read checking every base: per-base quality, average quality, and any <strong>N</strong>. Reads failing here are <strong>never aligned</strong>. Our example read is clean — it passes through."
}, async function () {
  label(500, 40, "scanning every base of the read", { anchor: "middle", size: 13, weight: 700 });

  const ry = 200;
  // the read as a row of nucleotide-tinted cells
  const span = [AMPLICON.fwd[1] + 1, AMPLICON.rev[0] - 1]; // between the read's primers
  const cellW = (bp2x(span[1]) - bp2x(span[0])) / (span[1] - span[0]);
  const readG = E("g", {}, canvas);
  for (let bp = span[0]; bp <= span[1]; bp++) {
    const nuc = "ACGT"[(bp * 7 + 3) % 4]; // deterministic fake sequence for tinting
    E("rect", { x: bp2x(bp), y: ry - 9, width: cellW + 0.5, height: 18, fill: nucColor(nuc), opacity: 0.55 }, readG);
  }
  label(500, ry - 22, "read", { anchor: "middle", size: 12, weight: 700 });

  await sleep(200);

  // scanning cursor sweeps left -> right
  const cursor = E("g", {}, canvas);
  E("rect", { x: -2, y: ry - 14, width: 4, height: 28, fill: "#0f172a" }, cursor);
  const cLab = E("text", { x: 0, y: ry - 20, "text-anchor": "middle", "font-size": 10, "font-weight": 700, fill: "#0f172a" }, cursor);
  cLab.textContent = "Q";
  place(cursor, bp2x(span[0]), 0);
  await anim(cursor, { x: bp2x(span[1]) }, 2200);

  // flag one position as 'N' to show what a failure looks like (faint, ghost)
  const npos = Math.round((span[0] + span[1]) / 2);
  E("circle", { cx: bp2x(npos), cy: ry, r: 0, fill: "none", stroke: "#ef4444", "stroke-width": 2 }, canvas);
  // (example read is clean, so we don't actually mark anything red — just demonstrate the scan)

  await sleep(200);
  // verdict
  const v = E("g", {}, canvas); place(v, 500, 320); setOpacity(v, 0);
  E("rect", { x: -90, y: -16, width: 180, height: 32, rx: 16, fill: "#dcfce7", stroke: "#22c55e" }, v);
  const vt = E("text", { x: 0, y: 5, "text-anchor": "middle", "font-size": 13, "font-weight": 700, fill: "#15803d" }, v);
  vt.textContent = "passes →";
  await fadeIn(v, 400);
});


/* 3. Collapse duplicates ----------------------------------------- */
defineScene({
  tag: "collapse",
  title: "Collapse identical reads",
  body: "Identical read pairs are collapsed into a single unique read carrying a <strong>×count</strong> weight. The same sequence is aligned only once; counts weight every downstream statistic — that is how 20 raw reads become 11 unique reads here."
}, async function () {
  label(500, 44, "Duplicate reads stack up and collapse into one weighted read", { anchor: "middle", size: 12, weight: 600 });

  // left: three identical reads
  const lefts = [];
  for (let i = 0; i < 3; i++) {
    const g = E("g", {}, canvas); place(g, 140, 110 + i * 40); setOpacity(g, 0);
    E("rect", { x: -55, y: -9, width: 110, height: 18, rx: 9, fill: "#fce7f3", stroke: "#ec4899" }, g);
    const lt = E("text", { x: 0, y: 3, "text-anchor": "middle", "font-size": 10, "font-weight": 600 }, g);
    lt.textContent = "same seq";
    lefts.push(g);
    fadeIn(g, 400);
    await sleep(180);
  }
  await sleep(250);
  // collapse arrow
  arrow(225, 170, 330, 170, "#94a3b8");
  label(278, 160, "collapse", { anchor: "middle", size: 10, fill: PAL.muted });
  await sleep(150);
  // merged
  const merged = E("g", {}, canvas); place(merged, 470, 170); setOpacity(merged, 0);
  E("rect", { x: -65, y: -12, width: 130, height: 24, rx: 12, fill: "#fce7f3", stroke: "#ec4899", "stroke-width": 1.6 }, merged);
  const mt = E("text", { x: 0, y: 4, "text-anchor": "middle", "font-size": 11, "font-weight": 700 }, merged);
  mt.textContent = "unique read";
  // count badge
  const badge = E("g", {}, merged); place(badge, 58, -14);
  E("circle", { cx: 0, cy: 0, r: 12, fill: "#ec4899" }, badge);
  const bt = E("text", { x: 0, y: 4, "text-anchor": "middle", "font-size": 11, "font-weight": 700, fill: "#fff" }, badge);
  bt.textContent = "×3";
  await fadeIn(merged, 500);

  // counters
  const cRaw = Counter(760, 150, "raw reads", 20);
  const cUniq = Counter(760, 250, "unique reads", 20, "#ec4899");
  await sleep(300);
  await cUniq.set(11, 700);
  label(760, 310, "(barcode_1)", { anchor: "middle", size: 10, fill: PAL.muted });
});

/* 4. Primer assignment ------------------------------------------- */
defineScene({
  tag: "assign",
  title: "Primer assignment",
  body: "Each read is assigned by locating its <strong>primers</strong> (overlap alignment, ~half-length minimum overlap, up to <code>primer_mismatch=2</code>). With the default <code>fastqfiles=0.5</code>, <strong>one</strong> primer match is enough — so the read is <em>assigned</em>. A read matching neither primer is <em>unassigned</em>."
}, async function () {
  label(500, 38, "find the config primers inside the read", { anchor: "middle", size: 13, weight: 700 });

  // amplicon reference at top
  const ay = 120;
  const amp = drawAmplicon(ay, { window: false });
  setOpacity(amp.g, 0); await fadeIn(amp.g, 500);
  label(bp2x((AMPLICON.fwd[0] + AMPLICON.fwd[1]) / 2), ay - 16, "fwd primer", { anchor: "middle", size: 9.5, fill: "#047857", weight: 600 });
  label(bp2x((AMPLICON.rev[0] + AMPLICON.rev[1]) / 2), ay - 16, "rev primer", { anchor: "middle", size: 9.5, fill: "#047857", weight: 600 });

  // two reads below, each with green primer ends
  const rA = 250, rB = 330;
  [[rA, "read A"], [rB, "read B"]].forEach(rd => {
    E("rect", { x: AX0, y: rd[0] - 6, width: AX1 - AX0, height: 12, rx: 6, fill: PAL.soft, stroke: PAL.line }, canvas);
    E("rect", { x: bp2x(AMPLICON.fwd[0]), y: rd[0] - 6, width: bp2x(AMPLICON.fwd[1]) - bp2x(AMPLICON.fwd[0]), height: 12, rx: 3, fill: PAL.primer }, canvas);
    E("rect", { x: bp2x(AMPLICON.rev[0]), y: rd[0] - 6, width: bp2x(AMPLICON.rev[1]) - bp2x(AMPLICON.rev[0]) + 2, height: 12, rx: 3, fill: PAL.primer }, canvas);
    label(AX0 - 8, rd[0] + 4, rd[1], { anchor: "end", size: 11, weight: 700 });
  });
  await sleep(200);

  // forward-primer probe slides from the amplicon down onto read A's left end and snaps
  const fwdW = bp2x(AMPLICON.fwd[1]) - bp2x(AMPLICON.fwd[0]);
  const probeF = E("rect", { x: 0, y: -5, width: fwdW, height: 10, rx: 3, fill: PAL.primer, opacity: 0.65 }, canvas);
  place(probeF, bp2x(AMPLICON.fwd[0]), ay);
  await sleep(150);
  await anim(probeF, { y: rA }, 550);
  await flashTick(bp2x(AMPLICON.fwd[0]) + fwdW / 2, rA - 16);

  // reverse-primer probe onto read A's right end
  const revW = bp2x(AMPLICON.rev[1]) - bp2x(AMPLICON.rev[0]) + 2;
  const probeR = E("rect", { x: 0, y: -5, width: revW, height: 10, rx: 3, fill: PAL.primer, opacity: 0.65 }, canvas);
  place(probeR, bp2x(AMPLICON.rev[0]), ay);
  await sleep(150);
  await anim(probeR, { y: rA }, 550);
  await flashTick(bp2x(AMPLICON.rev[0]) + revW / 2, rA - 16);

  await sleep(150);
  await badge(500, rA + 36, "assigned", "#dcfce7", "#15803d");

  // ghost read with neither primer matching -> unassigned
  const gy = 440;
  E("rect", { x: AX0, y: gy - 6, width: AX1 - AX0, height: 12, rx: 6, fill: "#f1f5f9", stroke: "#cbd5e1", "stroke-dasharray": "3 3" }, canvas);
  label(AX0 - 8, gy + 4, "?", { anchor: "end", size: 11, weight: 700, fill: PAL.muted });
  await sleep(150);
  await flashCross(bp2x((AMPLICON.fwd[0] + AMPLICON.fwd[1]) / 2), gy);
  await flashCross(bp2x((AMPLICON.rev[0] + AMPLICON.rev[1]) / 2), gy);
  await badge(500, gy + 34, "unassigned", "#fee2e2", "#b91c1c");

  label(500, 520, "fastqfiles = 0.5 → one primer match is enough", { anchor: "middle", size: 11, fill: PAL.muted });

  // ---- tiny local helpers ----------------------------------------
  async function flashTick(x, y) {
    const g = E("g", {}, canvas); place(g, x, y); setOpacity(g, 0);
    E("path", { d: "M-6,0 L-2,5 L6,-5", fill: "none", stroke: "#22c55e", "stroke-width": 3, "stroke-linecap": "round", "stroke-linejoin": "round" }, g);
    await fadeIn(g, 250);
  }
  async function flashCross(x, y) {
    const g = E("g", {}, canvas); place(g, x, y); setOpacity(g, 0);
    E("path", { d: "M-5,-5 L5,5 M5,-5 L-5,5", stroke: "#ef4444", "stroke-width": 3, "stroke-linecap": "round" }, g);
    await fadeIn(g, 250);
  }
  async function badge(x, y, txt, fill, stroke) {
    const g = E("g", {}, canvas); place(g, x, y); setOpacity(g, 0);
    E("rect", { x: -70, y: -14, width: 140, height: 28, rx: 14, fill: fill, stroke: stroke }, g);
    const t = E("text", { x: 0, y: 5, "text-anchor": "middle", "font-size": 13, "font-weight": 700, fill: stroke }, g);
    t.textContent = txt;
    await fadeIn(g, 400);
  }
});

/* 5. Global alignment ------------------------------------------- */
defineScene({
  tag: "align",
  title: "Global alignment",
  body: "The read and amplicon are <strong>anchored</strong> by their shared primers, then the middle is aligned with Needleman–Wunsch (<code>gap_opening=25</code>, <code>gap_extension=0</code>). A clean read matches end-to-end; an edited read opens a gap — the deletion. The pairwise views below are generated by R, with the same aligner amplican uses."
}, async function () {
  canvas.setAttribute("viewBox", "0 0 1000 300");
  label(500, 38, "primers anchor the alignment", { anchor: "middle", size: 13, weight: 700 });
  const ay = 110;
  const amp = drawAmplicon(ay);
  setOpacity(amp.g, 0); await fadeIn(amp.g, 500);

  // read starts offset to the right, then slides so its primers align under the amplicon's
  const ry = 220;
  const readG = E("g", {}, canvas);
  E("rect", { x: AX0, y: ry - 6, width: AX1 - AX0, height: 12, rx: 6, fill: PAL.soft, stroke: PAL.line }, readG);
  E("rect", { x: bp2x(AMPLICON.fwd[0]), y: ry - 6, width: bp2x(AMPLICON.fwd[1]) - bp2x(AMPLICON.fwd[0]), height: 12, rx: 3, fill: PAL.primer }, readG);
  E("rect", { x: bp2x(AMPLICON.rev[0]), y: ry - 6, width: bp2x(AMPLICON.rev[1]) - bp2x(AMPLICON.rev[0]) + 2, height: 12, rx: 3, fill: PAL.primer }, readG);
  place(readG, 80, 0); setOpacity(readG, 0);
  await anim(readG, { x: 0, opacity: 1 }, 650);

  // green dashed connectors between the matched primers
  [AMPLICON.fwd, AMPLICON.rev].forEach(p => {
    const x = bp2x((p[0] + p[1]) / 2);
    E("line", { x1: x, y1: ay + 8, x2: x, y2: ry - 7, stroke: "#86efac", "stroke-width": 1.6, "stroke-dasharray": "3 3" }, canvas);
  });
  await sleep(400);

  label(500, ry + 40, "two outcomes ↓", { anchor: "middle", size: 12, weight: 700, fill: PAL.muted });
  // embed the R-generated pairwise alignments in the figure pane
  showFigure("amplican-figures/alignment_match.svg", "matching read — no edits");
  await sleep(700);
  showFigure("amplican-figures/alignment_deletion.svg", "edited read — a deletion opens a gap (red)");
});

/* 6. HDR -------------------------------------------------------- */
defineScene({
  tag: "hdr",
  title: "HDR detection",
  body: "When a <strong>Donor</strong> is supplied, each read is checked for the donor-specific changes. Here the donor carries G→A; the read agrees at that position → <code>readType = HDR</code>."
}, async function () {
  label(500, 40, "donor vs amplicon vs read", { anchor: "middle", size: 13, weight: 700 });
  const startX = 270, trackLen = 460, yA = 150, yD = 230, yR = 310;
  const px = i => startX + (i / 30) * trackLen;
  [["amplicon", yA, "#e2e8f0"], ["donor", yD, "#bae6fd"], ["read", yR, "#fce7f3"]].forEach(t => {
    E("rect", { x: startX, y: t[1] - 9, width: trackLen, height: 18, rx: 5, fill: t[2], stroke: PAL.line }, canvas);
    label(startX - 96, t[1] + 4, t[0], { size: 11, weight: 700 });
  });
  const hx = px(15);
  E("rect", { x: hx - 11, y: yA - 16, width: 22, height: yR - yA + 30, fill: "#fde68a", opacity: 0.5 }, canvas);
  await sleep(200);
  [["G", yA, PAL.G], ["A", yD, PAL.A], ["A", yR, PAL.A]].forEach(l => {
    const c = E("circle", { cx: hx, cy: l[1], r: 0, fill: l[2] }, canvas);
    anim(c, { r: 9 }, 400);
    const t = E("text", { x: hx, y: l[1] + 4, "text-anchor": "middle", "font-size": 11, "font-weight": 700, fill: "#fff" }, canvas);
    t.textContent = l[1] === yA ? "G" : "A";
    setOpacity(t, 0); fadeIn(t, 400);
  });
  await sleep(500);
  arrow(hx, yA + 12, hx, yD - 12, "#0ea5e9");
  arrow(hx, yD + 12, hx, yR - 12, "#0ea5e9");
  label(hx + 80, (yA + yD) / 2, "donor change", { size: 10, fill: "#0369a1", weight: 600 });
  label(hx + 80, (yD + yR) / 2, "read = donor", { size: 10, fill: "#0369a1", weight: 600 });

  const v = E("g", {}, canvas); place(v, 500, 410); setOpacity(v, 0);
  E("rect", { x: -115, y: -18, width: 230, height: 36, rx: 18, fill: "#e0f2fe", stroke: PAL.hdr }, v);
  const vt = E("text", { x: 0, y: 6, "text-anchor": "middle", "font-size": 14, "font-weight": 700, fill: "#0369a1" }, v);
  vt.textContent = "→ readType = HDR";
  await fadeIn(v, 500);
});

/* 7. Event extraction ------------------------------------------- */
defineScene({
  tag: "extract",
  title: "Extract events",
  body: "The aligned strings are scanned to emit discrete <strong>events</strong> — deletions, insertions, mismatches — each with its position, type and strand. These are what every downstream plot counts."
}, async function () {
  canvas.setAttribute("viewBox", "0 0 1000 330");
  label(500, 38, "one alignment → a set of events", { anchor: "middle", size: 13, weight: 700 });
  const amp = drawAmplicon(110);
  setOpacity(amp.g, 0); await fadeIn(amp.g, 500);
  const rd = detailedRead(ID1_READS[0], 190, { fill: "#cbd5e1" });
  setOpacity(rd, 0); await fadeIn(rd, 500);
  await sleep(250);

  const evs = [
    { mid: (42 + 107) / 2, color: PAL.del, text: "deletion" },
    { mid: (108 + 127) / 2, color: PAL.ins, text: "insertion" },
    { mid: 31, color: nucColor("G"), text: "mismatch" },
    { mid: 163, color: nucColor("A"), text: "mismatch" }
  ];
  const chipY = 300;
  for (const ev of evs) {
    const cx = bp2x(ev.mid);
    const dot = E("circle", { cx: cx, cy: 190, r: 0, fill: ev.color }, canvas);
    anim(dot, { r: 5 }, 300);
    E("line", { x1: cx, y1: 196, x2: cx, y2: chipY - 14, stroke: ev.color, "stroke-width": 1, "stroke-dasharray": "2 2" }, canvas);
    const g = E("g", {}, canvas); place(g, cx, chipY); setOpacity(g, 0);
    E("rect", { x: -46, y: -13, width: 92, height: 26, rx: 13, fill: "#fff", stroke: ev.color, "stroke-width": 1.4 }, g);
    const t = E("text", { x: 0, y: 4, "text-anchor": "middle", "font-size": 11, "font-weight": 700, fill: ev.color }, g);
    t.textContent = ev.text;
    await fadeIn(g, 300);
    await sleep(180);
  }

  const efig = showFigure("amplican-figures/event_level_resolution.svg",
    "event-level resolution — deletion / insertion / mismatch");
  efig.style.width = "50%";
  efig.style.marginLeft = "auto";
  efig.style.marginRight = "auto";
  efig.style.marginTop = "-12px";
});


/* 8. Artifact filtering ----------------------------------------- */
defineScene({
  tag: "filter",
  title: "Filter artifacts",
  body: "Three filters remove non-biological events: <strong>EOP</strong> drops events overlapping primers; <strong>PRIMER_DIMER</strong> drops reads with a deletion wider than <code>amplicon − primers − 30</code>; <strong>findLQR</strong> clusters score-vs-events and drops the off-target blob (needs ≥1000 events)."
}, async function () {
  label(500, 32, "three artifact filters", { anchor: "middle", size: 13, weight: 700 });

  // ---- panel a: EOP ---------------------------------------------
  await panelTitle(60, "a)  events overlapping primers (EOP)");
  const a1 = drawAmplicon(96, { window: false });
  E("rect", { x: bp2x(160), y: 93, width: bp2x(199) - bp2x(160), height: 6, rx: 1, fill: PAL.del, opacity: 0.85 }, canvas);
  label(bp2x(180), 86, "deletion reaches the primer", { anchor: "middle", size: 9.5, fill: PAL.del, weight: 700 });
  await removedBadge(bp2x(180), 124, "event removed");

  // ---- panel b: primer dimer -----------------------------------
  await panelTitle(168, "b)  primer dimer  (deletion > amplicon − primers − 30)");
  drawAmplicon(204, { window: false });
  E("rect", { x: bp2x(25), y: 201, width: bp2x(175) - bp2x(25), height: 6, rx: 1, fill: PAL.del, opacity: 0.85 }, canvas);
  label(bp2x(100), 194, "150 bp deletion", { anchor: "middle", size: 9.5, fill: PAL.del, weight: 700 });
  await removedBadge(bp2x(100), 232, "read removed");

  // ---- panel c: off-target / low-quality (findLQR) -------------
  await panelTitle(276, "c)  off-target / low-quality  (findLQR, CLARA clustering)");
  const sx = 300, sy = 306, sw = 400, sh = 150;
  E("rect", { x: sx, y: sy, width: sw, height: sh, fill: "#fff", stroke: PAL.line }, canvas);
  label(sx + sw, sy + sh + 18, "alignment score  →", { anchor: "end", size: 9.5, fill: PAL.muted });
  label(sx + 10, sy - 6, "events", { size: 9.5, fill: PAL.muted, weight: 700 });
  for (let i = 0; i < 30; i++) {
    const good = i % 4 !== 0;
    const cx2 = sx + (good ? 90 + Math.random() * 300 : 15 + Math.random() * 55);
    const cy2 = sy + sh - (good ? 15 + Math.random() * 55 : 85 + Math.random() * 50);
    const d = E("circle", { cx: cx2, cy: cy2, r: 3, fill: good ? "#86efac" : "#fca5a5", opacity: 0 }, canvas);
    later(() => anim(d, { opacity: 1 }, 200), i * 25);
  }
  await sleep(700);
  E("circle", { cx: sx + 40, cy: sy + sh - 95, r: 38, fill: "none", stroke: "#ef4444", "stroke-width": 1.5, "stroke-dasharray": "4 3" }, canvas);
  label(sx + 150, sy + 18, "off-target cluster → dropped", { size: 10, fill: PAL.del, weight: 700 });
  label(sx + sw - 10, sy + sh - 8, "activates only with ≥ 1000 events", { anchor: "end", size: 9, fill: PAL.muted });

  // ---- local helpers -------------------------------------------
  async function panelTitle(y, txt) {
    const t = label(120, y, txt, { size: 11.5, weight: 700, fill: "#334155" });
    setOpacity(t, 0); await fadeIn(t, 250);
  }
  async function removedBadge(x, y, txt) {
    const g = E("g", {}, canvas); place(g, x, y); setOpacity(g, 0);
    E("rect", { x: -58, y: -11, width: 116, height: 22, rx: 11, fill: "#fee2e2", stroke: "#ef4444" }, g);
    const t = E("text", { x: 0, y: 4, "text-anchor": "middle", "font-size": 10.5, "font-weight": 700, fill: "#b91c1c" }, g);
    t.textContent = "✕ " + txt;
    await fadeIn(g, 350);
  }
});

/* 9. Cut-site overlap ------------------------------------------- */
defineScene({
  tag: "overlap",
  title: "Cut-site overlap",
  body: "Events overlapping the cut-site window (UPPERCASE ± <code>cut_buffer=5</code>, here <strong>61–88</strong>) are flagged as CRISPR edits. The deletion 42–107 <strong>overlaps</strong> (pink, counts); the insertion 108–127 sits <strong>outside</strong> (blue, does not count). That is why ID_1 has Reads_Del=6 but Reads_In=0."
}, async function () {
  label(500, 40, "Inside the window → counted. Outside → ignored.", { anchor: "middle", size: 12, weight: 600 });
  // Create the window first so the amplicon and its labels stay above it.
  const wBand = E("rect", { x: bp2x(61), y: 90, width: bp2x(88) - bp2x(61), height: 80, fill: PAL.cutwin, opacity: 0.6 }, canvas);
  setOpacity(wBand, 0);

  const amp = drawAmplicon(120);
  annotateAmplicon(amp, true);

  await fadeIn(wBand, 400);
  label(bp2x(74.5), 86, "window 61–88", { anchor: "middle", size: 10, weight: 700, fill: "#b45309" });

  const rd = detailedRead(ID1_READS[0], 180, { fill: "#cbd5e1" });
  setOpacity(rd, 0); await fadeIn(rd, 500);

  await sleep(300);
  // deletion (overlaps) -> pink
  const d = E("rect", { x: bp2x(42), y: 167, width: 0, height: 6, fill: PAL.pink }, canvas);
  await anim(d, { width: bp2x(107) - bp2x(42) }, 600);
  label(bp2x(74), 162, "deletion overlaps → counts", { anchor: "middle", size: 10, fill: PAL.pink, weight: 700 });

  // insertion (outside) -> blue/grey
  const px = bp2x((108 + 127) / 2);
  label(px + 12, 196, "insertion 108–127 → outside, ignored", { size: 10, fill: PAL.ins, weight: 700 });

  const cDel = Counter(250, 460, "Reads_Del", 0, PAL.pink);
  const cIn = Counter(500, 460, "Reads_In", 0, PAL.ins);
  const cEd = Counter(750, 460, "Reads_Edited", 0, "#6d28d9");
  await sleep(300);
  await cDel.set(6, 600);
  await cIn.set(0, 400);
  await cEd.set(6, 600);
});

/* 10. Relative coordinates -------------------------------------- */
defineScene({
  tag: "shift",
  title: "Shift to relative coords",
  body: "Positions are re-zeroed at the cut site: upstream negative, downstream positive. This lets <strong>metaplots</strong> overlay different amplicons. For <code>Direction=1</code> amplicons the strand is flipped so everything is reported on the forward strand."
}, async function () {
  label(500, 40, "absolute → cut-site-relative", { anchor: "middle", size: 12, weight: 600 });
  const amp = drawAmplicon(140);

  // absolute ruler
  label(60, 100, "absolute (bp)", { size: 11, weight: 700 });
  [1, 30, 66, 100, 140, 178].forEach(bp => {
    E("line", { x1: bp2x(bp), y1: 110, x2: bp2x(bp), y2: 118, stroke: PAL.muted }, canvas);
    label(bp2x(bp), 108, bp, { anchor: "middle", size: 9 });
  });

  await sleep(400);
  // cut site zero marker
  E("line", { x1: bp2x(66), y1: 90, x2: bp2x(66), y2: 300, stroke: PAL.cut, "stroke-dasharray": "3 3" }, canvas);

  // relative ruler animates numbers
  label(60, 230, "relative (cut = 0)", { size: 11, weight: 700, fill: "#b45309" });
  const rels = [[1, -65], [30, -36], [66, 0], [100, 34], [140, 74], [178, 112]];
  const numEls = [];
  rels.forEach(r => {
    const t = label(bp2x(r[0]), 218, r[0], { anchor: "middle", size: 9, fill: PAL.muted });
    numEls.push({ el: t, from: r[0], to: r[1] });
    E("line", { x1: bp2x(r[0]), y1: 228, x2: bp2x(r[0]), y2: 236, stroke: PAL.cut }, canvas);
  });
  // tween the numbers
  const t0 = performance.now();
  await new Promise(r => {
    function f(now) {
      const t = Math.min(1, (now - t0) / 900);
      const e = easeInOut(t);
      numEls.forEach(n => n.el.textContent = Math.round(n.from + (n.to - n.from) * e));
      if (t < 1) requestAnimationFrame(f); else r();
    }
    requestAnimationFrame(f);
  });

  await sleep(200);
  label(bp2x(66) + 6, 270, "0 = first UPPERCASE base", { size: 10, fill: "#b45309", weight: 700 });
});

/* 11. Normalization --------------------------------------------- */
defineScene({
  tag: "normalize",
  title: "Normalize vs controls",
  body: "Control samples define background. Events in treatments that <strong>exactly match</strong> a control event (same pos/type/base, above <code>min_freq=1%</code>) are removed; control events are kept. Handles SNPs and index-hopping."
}, async function () {
  canvasWrap.style.display = "none";
  showFigure("amplican-figures/amplicanNormalize.svg",
    "amplicanNormalize — control background is subtracted from treatments");
});

/* 12. Summarization --------------------------------------------- */
defineScene({
  tag: "summarize",
  title: "Summarize per experiment",
  body: "Per-read outcomes roll up to per-experiment tallies (weighted by count). Only <strong>consensus</strong> events <strong>overlapping the cut site</strong> count. Frameshift = net indel (cut-overlapping) not divisible by 3 → read 2's 74 bp deletion is a frameshift (×2)."
}, async function () {
  label(500, 36, "Reads → counters (ID_1)", { anchor: "middle", size: 12, weight: 600 });

  // left: stack of reads with counts
  const reads = [
    { id: "read 3 ×1", note: "primer dimer", fate: "dimer" },
    { id: "read 1 ×3", note: "del 66 (overlap)", fate: "del" },
    { id: "read 2 ×2", note: "del 74 → frameshift", fate: "fs" },
    { id: "read 4 ×1", note: "del 84 (overlap)", fate: "del" }
  ];
  const lx = 60;
  reads.forEach((r, i) => {
    const y = 80 + i * 40;
    const g = E("g", {}, canvas); place(g, lx, y); setOpacity(g, 0);
    E("rect", { x: 0, y: -12, width: 230, height: 28, rx: 8, fill: "#fff", stroke: PAL.line }, g);
    const t1 = E("text", { x: 10, y: 6, "font-size": 11, "font-weight": 700 }, g);
    t1.textContent = r.id;
    const c = r.fate === "dimer" ? "#8b5cf6" : r.fate === "fs" ? "#ef4444" : PAL.pink;
    const t2 = E("text", { x: 120, y: 6, "font-size": 10, fill: c }, g);
    t2.textContent = r.note;
    later(() => fadeIn(g, 400), i * 250);
  });

  // right: counters grid tally up
  const grid = [
    ["Reads", 7], ["PRIMER_DIMER", 1], ["Low_Score", 0],
    ["Reads_Filtered", 6], ["Reads_Del", 6], ["Reads_In", 0],
    ["Reads_Edited", 6], ["Frameshifted", 2], ["HDR", 0]
  ];
  const cs = [];
  const gx = 380, gy = 80;
  grid.forEach((cell, i) => {
    const r = Math.floor(i / 3), c = i % 3;
    const x = gx + c * 200, y = gy + r * 90;
    const cnt = Counter(x, y, cell[0], 0,
      ["Reads_Del", "Frameshifted"].includes(cell[0]) ? PAL.pink :
      cell[0] === "HDR" ? PAL.hdr : cell[0] === "PRIMER_DIMER" ? "#8b5cf6" : PAL.ink);
    cs.push({ cnt, v: cell[1] });
  });

  await sleep(600);
  for (const c of cs) {
    await c.cnt.set(c.v, 500);
    await sleep(120);
  }

  // frameshift math footnote
  label(500, 410, "frameshift math:", { anchor: "middle", size: 11, weight: 700 });
  label(500, 432, "read 2: deletion 74 bp overlaps cut  →  net indel −74  →  −74 mod 3 ≠ 0  →  frameshift (×2)", { anchor: "middle", size: 10.5, fill: PAL.del });

  // common rate equations (formulas only, no numbers)
  label(500, 466, "common rates:", { anchor: "middle", size: 11, weight: 700 });
  label(500, 488, "edit rate % = round(Reads_Edited × 100 / Reads_Filtered, 2)", { anchor: "middle", size: 10.5, fill: PAL.pink });
  label(500, 510, "HDR % = round(HDR × 100 / Reads_Filtered, 2)", { anchor: "middle", size: 10.5, fill: PAL.hdr });
});

/* 13. Reporting ------------------------------------------------- */
defineScene({
  tag: "reports",
  title: "Reports",
  body: "Six HTML reports are knitted. <strong>Index</strong> is the landing page — data-quality checks and links to everything else. <strong>Barcode</strong> verifies read assignment per pool. <strong>Guide</strong>, <strong>Group</strong> and <strong>Amplicon</strong> are <em>metaplots</em> — they aggregate events on the shared cut-relative axis, stratified by that column. <strong>ID</strong> is the per-sample ground truth."
}, async function () {
  label(500, 28, "Six HTML reports, four roles", { anchor: "middle", size: 12, weight: 600 });

  // ---- role motifs (no arches) ---------------------------------
  function motifIndex(g, w) {
    const xs = [w / 2 - 20, w / 2, w / 2 + 20];
    const cols = ["#22c55e", "#22c55e", "#f59e0b"];
    xs.forEach((x, i) => E("circle", { cx: x, cy: 52, r: 4, fill: cols[i] }, g));
  }
  function motifBarcode(g, w) {
    E("rect", { x: w / 2 - 34, y: 44, width: 68, height: 16, rx: 8, fill: "#e0f2fe", stroke: "#0ea5e9" }, g);
    const t = E("text", { x: w / 2, y: 56, "text-anchor": "middle", "font-size": 10, "font-weight": 600, fill: "#0369a1" }, g);
    t.textContent = "reads  \u2713";
  }
  function motifAgg(g, w) {
    const baseY = 92;
    const xs = [w / 2 - 24, w / 2, w / 2 + 24];
    const hs = [20, 34, 26];
    const cols = [PAL.pink, "#60a5fa", PAL.pink];
    xs.forEach((x, i) => E("rect", { x: x - 5, y: baseY - hs[i], width: 10, height: hs[i], rx: 2, fill: cols[i], opacity: 0.8 }, g));
  }
  function motifID(g, w) {
    E("rect", { x: w / 2 - 46, y: 46, width: 92, height: 16, rx: 8, fill: "#faf5ff", stroke: "#a855f7" }, g);
    const t = E("text", { x: w / 2, y: 58, "text-anchor": "middle", "font-size": 10, "font-weight": 600, fill: "#6b21a8" }, g);
    t.textContent = "1 sample";
    E("circle", { cx: w / 2 + 30, cy: 54, r: 4, fill: PAL.A, stroke: "#fff", "stroke-width": 1 }, g);
  }
  function rNode(cx, y, w, h, title, sub, color, motif) {
    const x = cx - w / 2;
    const g = E("g", {}, canvas); place(g, x, y); setOpacity(g, 0);
    E("rect", { x: 0, y: 0, width: w, height: h, rx: 10, fill: "#fff", stroke: color, "stroke-width": 1.6 }, g);
    E("rect", { x: 0, y: 0, width: w, height: 24, rx: 10, fill: color }, g);
    E("rect", { x: 0, y: 14, width: w, height: 10, fill: color }, g);
    const t = E("text", { x: w / 2, y: 16, "text-anchor": "middle", "font-size": 12, "font-weight": 700, fill: "#fff" }, g);
    t.textContent = title;
    const s = E("text", { x: w / 2, y: 42, "text-anchor": "middle", "font-size": 10.5, fill: PAL.muted }, g);
    s.textContent = sub;
    if (motif) motif(g, w);
    return g;
  }

  // Index — hub
  const indexN = rNode(500, 44, 260, 64, "Index", "landing \u2014 data quality & links to all", "#ec4899", motifIndex);
  await fadeIn(indexN, 400); await sleep(160);
  arrow(500, 108, 500, 124);

  // Barcode — QC gate
  const barN = rNode(500, 124, 260, 64, "Barcode", "read-assignment checks (per pool)", "#0ea5e9", motifBarcode);
  await fadeIn(barN, 400); await sleep(160);

  // bus Index/Barcode level -> metaplot trio
  E("line", { x1: 500, y1: 188, x2: 500, y2: 208, stroke: PAL.muted, "stroke-width": 1.4 }, canvas);
  E("line", { x1: 200, y1: 208, x2: 800, y2: 208, stroke: PAL.muted, "stroke-width": 1.4 }, canvas);
  [200, 500, 800].forEach(x => E("line", { x1: x, y1: 208, x2: x, y2: 230, stroke: PAL.muted, "stroke-width": 1.4 }, canvas));
  // metaplots bracket + rotated label
  E("path", { d: "M70,230 L54,230 L54,334 L70,334", fill: "none", stroke: PAL.muted, "stroke-width": 1.4 }, canvas);
  const mlbl = E("text", { x: 38, y: 282, "text-anchor": "middle", "font-size": 10, "font-weight": 700, fill: PAL.muted, transform: "rotate(-90 38 282)" }, canvas);
  mlbl.textContent = "metaplots";

  // trio: Guide / Group / Amplicon
  const guideN = rNode(200, 230, 240, 104, "Guide", "per guide", "#f59e0b", motifAgg);
  const groupN = rNode(500, 230, 240, 104, "Group", "per group", "#10b981", motifAgg);
  const ampN = rNode(800, 230, 240, 104, "Amplicon", "per amplicon", "#6d28d9", motifAgg);
  await fadeIn(guideN, 350); await sleep(120);
  await fadeIn(groupN, 350); await sleep(120);
  await fadeIn(ampN, 350); await sleep(200);

  // bus trio -> ID
  [200, 500, 800].forEach(x => E("line", { x1: x, y1: 334, x2: x, y2: 366, stroke: PAL.muted, "stroke-width": 1.4 }, canvas));
  E("line", { x1: 200, y1: 366, x2: 800, y2: 366, stroke: PAL.muted, "stroke-width": 1.4 }, canvas);
  arrow(500, 366, 500, 392);

  // ID — per-sample truth
  const idN = rNode(500, 392, 300, 80, "ID report", "per-experiment ground truth", "#a855f7", motifID);
  await fadeIn(idN, 400); await sleep(200);

  label(500, 500, "done \u2014 from FASTQ to quantification", { anchor: "middle", size: 12, weight: 700, fill: "#6d28d9" });
});

/* ------------------------------------------------------------------ */
/* Boot                                                                */
/* ------------------------------------------------------------------ */
buildChapters();
show(0);

})();

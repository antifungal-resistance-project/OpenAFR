"use client";

import {
  useMotionValueEvent,
  useReducedMotion,
  useScroll,
} from "motion/react";
import { useRef, useState } from "react";

/*
 * The signature scroll animation.
 *
 * A candidate molecule approaches CYP51's heme along the ray through the iron.
 * As the reader scrolls, its nitrogen closes on the Fe and the N–Fe distance
 * ticks from ~4.6 Å down to 2.63 Å — the iron-coordinating pose that OpenAFR's
 * validated criterion rewards, matching the VT-1161 redock (work/RESULTS_redock_VT1.md).
 *
 * Geometry is derived from a single scroll progress value so everything stays
 * in lockstep. Reduced-motion users get the resolved end state, no scroll pin.
 */

const FE = { x: 200, y: 200 };
const DIR = { x: 120, y: -120 }; // Fe + t*DIR traces the approach ray (t: 1 → far)
const T_START = 1;
const T_END = 0.32;
const A_FAR = 4.62;
const A_BOUND = 2.63;

function clamp01(v: number) {
  return Math.min(1, Math.max(0, v));
}

function frame(t: number) {
  const nx = FE.x + DIR.x * t;
  const ny = FE.y + DIR.y * t;
  const dist =
    A_BOUND + ((t - T_END) / (T_START - T_END)) * (A_FAR - A_BOUND);
  // away-from-iron unit direction, for placing the ring off the coordinating N
  const away = { x: 0.7071, y: -0.7071 };
  return { nx, ny, dist, away };
}

function hexPath(cx: number, cy: number, r: number) {
  const pts: string[] = [];
  for (let i = 0; i < 6; i++) {
    const a = (Math.PI / 3) * i - Math.PI / 6;
    pts.push(`${(cx + r * Math.cos(a)).toFixed(1)},${(cy + r * Math.sin(a)).toFixed(1)}`);
  }
  return `M${pts.join("L")}Z`;
}

const CAPTIONS = [
  {
    at: 0.05,
    title: "Meet the target",
    body: "Deep inside a fungal cell sits CYP51 — an enzyme the fungus needs to survive. At its heart is an iron atom, held in a ring called a heme.",
  },
  {
    at: 0.33,
    title: "How the drugs work",
    body: "Antifungal drugs slip into this pocket and reach a nitrogen atom toward the iron, jamming the enzyme. That's how the whole azole class kills fungus.",
  },
  {
    at: 0.6,
    title: "Why resistance is hard to beat",
    body: "Many molecules look like they should fit. What actually matters is a matter of ångströms: how close that nitrogen truly gets to the iron.",
  },
  {
    at: 0.85,
    title: "The line between hit and miss",
    body: "N–Fe ≈ 2.63 Å — an iron-coordinating pose, the hallmark of a drug that works. Measuring this is how our research tells real candidates from look-alikes.",
  },
];

function activeCaption(p: number) {
  let idx = 0;
  for (let i = 0; i < CAPTIONS.length; i++) if (p >= CAPTIONS[i].at) idx = i;
  return idx;
}

export function HemeApproach() {
  const ref = useRef<HTMLDivElement>(null);
  const reduce = useReducedMotion();
  const { scrollYProgress } = useScroll({
    target: ref,
    offset: ["start start", "end end"],
  });

  const [p, setP] = useState(reduce ? 1 : 0);
  useMotionValueEvent(scrollYProgress, "change", (v) => {
    if (!reduce) setP(v);
  });

  // Map scroll progress → approach parameter t (held at the ends).
  const approach = clamp01((p - 0.12) / (0.9 - 0.12));
  const t = T_START + (T_END - T_START) * approach;
  const { nx, ny, dist, away } = frame(t);
  const ringCx = nx + away.x * 30;
  const ringCy = ny + away.y * 30;
  const bondOpacity = clamp01((p - 0.6) / 0.25);
  const glow = 0.5 + 0.5 * Math.sin(p * Math.PI * 4);
  const activeIdx = activeCaption(p);

  return (
    <section ref={ref} className="relative h-[320vh]">
      <div className="sticky top-0 flex h-screen items-center overflow-hidden">
        <div className="glow-iron pointer-events-none absolute inset-0" />
        <div className="mx-auto grid w-full max-w-6xl grid-cols-1 items-center gap-8 px-6 md:grid-cols-2">
          {/* Narration */}
          <div className="order-2 md:order-1">
            <p className="font-mono text-xs uppercase tracking-[0.2em] text-iron">
              Follow the molecule
            </p>
            <div className="relative mt-4 min-h-[180px]">
              {CAPTIONS.map((c, i) => (
                <div
                  key={c.title}
                  className="absolute inset-0 transition-all duration-500"
                  style={{
                    opacity: i === activeIdx ? 1 : 0.18,
                    transform: `translateY(${(i - activeIdx) * 8}px)`,
                  }}
                  aria-hidden={i !== activeIdx}
                >
                  <h3 className="text-2xl font-semibold tracking-tight text-fg sm:text-3xl">
                    {c.title}
                  </h3>
                  <p className="mt-3 max-w-md text-base leading-relaxed text-muted">
                    {c.body}
                  </p>
                </div>
              ))}
            </div>
            <div className="mt-8 flex items-baseline gap-3">
              <span
                className="font-mono text-5xl font-semibold tabular-nums tracking-tight sm:text-6xl"
                style={{ color: dist <= 3.0 ? "var(--color-pass)" : "var(--color-nitrogen)" }}
              >
                {dist.toFixed(2)}
              </span>
              <span className="font-mono text-lg text-muted">Å&nbsp;N–Fe</span>
            </div>
          </div>

          {/* Molecular view */}
          <div className="order-1 md:order-2">
            <svg
              viewBox="0 0 400 400"
              className="mx-auto w-full max-w-md"
              role="img"
              aria-label="A candidate molecule's nitrogen approaching the heme iron of CYP51"
            >
              {/* iron glow */}
              <circle
                cx={FE.x}
                cy={FE.y}
                r={46}
                fill="var(--color-iron)"
                opacity={0.1 + glow * 0.12}
              />
              {/* porphyrin macrocycle */}
              <circle
                cx={FE.x}
                cy={FE.y}
                r={70}
                fill="none"
                stroke="var(--color-border)"
                strokeWidth={2}
              />
              {[0, 90, 180, 270].map((deg) => {
                const a = (deg * Math.PI) / 180;
                const px = FE.x + 26 * Math.cos(a);
                const py = FE.y + 26 * Math.sin(a);
                const ox = FE.x + 70 * Math.cos(a);
                const oy = FE.y + 70 * Math.sin(a);
                return (
                  <g key={deg}>
                    <line
                      x1={FE.x}
                      y1={FE.y}
                      x2={px}
                      y2={py}
                      stroke="var(--color-faint)"
                      strokeWidth={1.5}
                    />
                    <line
                      x1={px}
                      y1={py}
                      x2={ox}
                      y2={oy}
                      stroke="var(--color-faint)"
                      strokeWidth={1.5}
                    />
                    <circle cx={px} cy={py} r={5} fill="var(--color-nitrogen)" opacity={0.55} />
                  </g>
                );
              })}
              {/* iron */}
              <circle cx={FE.x} cy={FE.y} r={12} fill="var(--color-iron)" />
              <text
                x={FE.x}
                y={FE.y + 4}
                textAnchor="middle"
                className="font-mono"
                fontSize={11}
                fill="#0a0e1a"
                fontWeight={700}
              >
                Fe
              </text>

              {/* coordinating bond, forms as the nitrogen closes in */}
              <line
                x1={FE.x}
                y1={FE.y}
                x2={nx}
                y2={ny}
                stroke="var(--color-pass)"
                strokeWidth={2}
                strokeDasharray="4 4"
                opacity={bondOpacity}
              />

              {/* candidate molecule */}
              <g>
                <path
                  d={hexPath(ringCx, ringCy, 24)}
                  fill="none"
                  stroke="var(--color-nitrogen)"
                  strokeWidth={2.5}
                  opacity={0.9}
                />
                <path
                  d={hexPath(ringCx, ringCy, 15)}
                  fill="none"
                  stroke="var(--color-nitrogen)"
                  strokeWidth={1.2}
                  opacity={0.5}
                />
                <circle cx={nx} cy={ny} r={9} fill="var(--color-nitrogen)" />
                <text
                  x={nx}
                  y={ny + 3.5}
                  textAnchor="middle"
                  className="font-mono"
                  fontSize={10}
                  fill="#05070d"
                  fontWeight={700}
                >
                  N
                </text>
              </g>
            </svg>
          </div>
        </div>

        {/* scroll hint, fades once the reader engages */}
        <div
          className="pointer-events-none absolute inset-x-0 bottom-6 flex justify-center transition-opacity duration-500"
          style={{ opacity: clamp01(1 - p * 6) }}
        >
          <span className="font-mono text-xs uppercase tracking-[0.2em] text-faint">
            scroll ↓
          </span>
        </div>
      </div>
    </section>
  );
}

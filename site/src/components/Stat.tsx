"use client";

import {
  animate,
  useInView,
  useReducedMotion,
} from "motion/react";
import { useEffect, useRef, useState } from "react";

type Kind = "decimal" | "p" | "int";

function format(value: number, kind: Kind, suffix?: string) {
  let body: string;
  if (kind === "p") body = value.toFixed(4);
  else if (kind === "int") body = Math.round(value).toLocaleString();
  else body = value.toFixed(value < 10 ? 3 : 2);
  return (kind === "p" ? "p = " + body : body) + (suffix ?? "");
}

/**
 * Counts up to `value` when scrolled into view. Tone tints the number:
 * "pass" = validated (emerald), "fail" = worse-than-random (rose).
 */
export function Stat({
  value,
  label,
  kind = "decimal",
  suffix,
  tone = "pass",
}: {
  value: number;
  label: string;
  kind?: Kind;
  suffix?: string;
  tone?: "pass" | "fail" | "neutral";
}) {
  const ref = useRef<HTMLDivElement>(null);
  const inView = useInView(ref, { once: true, margin: "-15% 0px" });
  const reduce = useReducedMotion();
  const [display, setDisplay] = useState(() =>
    reduce ? format(value, kind, suffix) : format(0, kind, suffix),
  );

  useEffect(() => {
    if (!inView || reduce) {
      if (reduce) setDisplay(format(value, kind, suffix));
      return;
    }
    const controls = animate(0, value, {
      duration: 1.4,
      ease: [0.16, 1, 0.3, 1],
      onUpdate: (v) => setDisplay(format(v, kind, suffix)),
    });
    return () => controls.stop();
  }, [inView, reduce, value, kind, suffix]);

  const toneClass =
    tone === "pass"
      ? "text-pass"
      : tone === "fail"
        ? "text-fail"
        : "text-iron";

  return (
    <div
      ref={ref}
      className="rounded-xl border border-border bg-surface/60 p-6 backdrop-blur"
    >
      <div
        className={`font-mono text-4xl font-semibold tabular-nums tracking-tight sm:text-5xl ${toneClass}`}
      >
        {display}
      </div>
      <div className="mt-3 text-sm leading-snug text-muted">{label}</div>
    </div>
  );
}

"use client";

import { useEffect, useState } from "react";
import { org } from "@/content/site";

const links = [
  { href: "#mission", label: "Mission" },
  { href: "#approach", label: "What we do" },
  { href: "#science", label: "The science" },
  { href: "#work", label: "Our work" },
];

export function Nav() {
  const [scrolled, setScrolled] = useState(false);

  useEffect(() => {
    const onScroll = () => setScrolled(window.scrollY > 24);
    onScroll();
    window.addEventListener("scroll", onScroll, { passive: true });
    return () => window.removeEventListener("scroll", onScroll);
  }, []);

  return (
    <header
      className={`fixed inset-x-0 top-0 z-50 transition-colors duration-300 ${
        scrolled
          ? "border-b border-border bg-bg/80 backdrop-blur"
          : "border-b border-transparent"
      }`}
    >
      <nav className="mx-auto flex max-w-6xl items-center justify-between px-6 py-4">
        <a href="#top" className="flex items-center gap-2.5">
          <span className="h-2.5 w-2.5 rounded-full bg-iron shadow-[0_0_12px_var(--color-iron)]" />
          <span className="font-mono text-sm font-semibold tracking-tight">
            {org.short}
          </span>
          <span className="hidden text-sm text-muted sm:inline">
            {org.name}
          </span>
        </a>
        <div className="hidden items-center gap-7 md:flex">
          {links.map((l) => (
            <a
              key={l.href}
              href={l.href}
              className="text-sm text-muted transition-colors hover:text-fg"
            >
              {l.label}
            </a>
          ))}
        </div>
        <a
          href="#contact"
          className="rounded-full border border-border px-4 py-1.5 text-sm text-fg transition-colors hover:border-iron hover:text-iron"
        >
          Get involved
        </a>
      </nav>
    </header>
  );
}

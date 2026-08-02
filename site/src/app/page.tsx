import { Nav } from "@/components/Nav";
import { Reveal } from "@/components/Reveal";
import { Stat } from "@/components/Stat";
import { HemeApproach } from "@/components/HemeApproach";
import {
  approach,
  cta,
  hero,
  org,
  problem,
  project,
  projectSection,
  roadmap,
  science,
  values,
} from "@/content/site";

function Eyebrow({ children }: { children: React.ReactNode }) {
  return (
    <p className="font-mono text-xs uppercase tracking-[0.2em] text-iron">
      {children}
    </p>
  );
}

const stateStyles: Record<string, { dot: string; label: string; text: string }> = {
  done: { dot: "bg-pass", label: "Done", text: "text-pass" },
  next: { dot: "bg-iron", label: "Underway", text: "text-iron" },
  blocked: { dot: "bg-fail", label: "Seeking partners", text: "text-fail" },
};

export default function Home() {
  return (
    <>
      <Nav />
      <main id="top">
        {/* Hero — the mission */}
        <section className="relative flex min-h-screen items-center overflow-hidden">
          <div className="bg-grid pointer-events-none absolute inset-0 opacity-60" />
          <div className="glow-iron pointer-events-none absolute inset-0" />
          <div className="absolute inset-x-0 bottom-0 h-40 bg-gradient-to-b from-transparent to-bg" />
          <div className="relative mx-auto max-w-4xl px-6 pt-28 pb-20">
            <Reveal>
              <Eyebrow>{hero.eyebrow}</Eyebrow>
            </Reveal>
            <Reveal delay={0.08}>
              <h1 className="mt-6 text-balance text-4xl font-semibold leading-[1.08] tracking-tight sm:text-6xl">
                {hero.title}
              </h1>
            </Reveal>
            <Reveal delay={0.16}>
              <p className="mt-6 max-w-2xl text-pretty text-lg leading-relaxed text-muted">
                {hero.sub}
              </p>
            </Reveal>
            <Reveal delay={0.24}>
              <div className="mt-9 flex flex-wrap items-center gap-4">
                <a
                  href="#mission"
                  className="rounded-full bg-iron px-6 py-3 text-sm font-semibold text-bg transition-transform hover:-translate-y-0.5"
                >
                  Our mission
                </a>
                <a
                  href="#contact"
                  className="rounded-full border border-border px-6 py-3 text-sm font-semibold text-fg transition-colors hover:border-iron hover:text-iron"
                >
                  Get involved
                </a>
              </div>
            </Reveal>
          </div>
        </section>

        {/* Problem — why fungi, why now */}
        <section id="mission" className="border-t border-border py-28">
          <div className="mx-auto max-w-6xl px-6">
            <Reveal>
              <Eyebrow>{problem.eyebrow}</Eyebrow>
            </Reveal>
            <Reveal delay={0.08}>
              <h2 className="mt-5 max-w-3xl text-3xl font-semibold tracking-tight sm:text-4xl">
                {problem.title}
              </h2>
            </Reveal>
            <div className="mt-14 grid gap-6 md:grid-cols-3">
              {problem.points.map((pt, i) => (
                <Reveal key={pt.k} delay={i * 0.1} as="div">
                  <div className="h-full rounded-2xl border border-border bg-surface/50 p-7">
                    <h3 className="text-lg font-semibold text-fg">{pt.k}</h3>
                    <p className="mt-3 text-[15px] leading-relaxed text-muted">
                      {pt.v}
                    </p>
                  </div>
                </Reveal>
              ))}
            </div>
          </div>
        </section>

        {/* Approach — what we do */}
        <section id="approach" className="border-t border-border py-28">
          <div className="mx-auto max-w-6xl px-6">
            <Reveal>
              <Eyebrow>{approach.eyebrow}</Eyebrow>
            </Reveal>
            <Reveal delay={0.08}>
              <h2 className="mt-5 max-w-3xl text-3xl font-semibold tracking-tight sm:text-4xl">
                {approach.title}
              </h2>
            </Reveal>
            <div className="mt-14 grid gap-6 md:grid-cols-3">
              {approach.pillars.map((pt, i) => (
                <Reveal key={pt.k} delay={i * 0.1} as="div">
                  <div className="h-full rounded-2xl border border-border bg-surface/50 p-7">
                    <div className="mb-5 font-mono text-sm text-iron">
                      {String(i + 1).padStart(2, "0")}
                    </div>
                    <h3 className="text-lg font-semibold text-fg">{pt.k}</h3>
                    <p className="mt-3 text-[15px] leading-relaxed text-muted">
                      {pt.v}
                    </p>
                  </div>
                </Reveal>
              ))}
            </div>
          </div>
        </section>

        {/* The science of resistance — scroll-driven animation */}
        <section id="science" className="border-t border-border">
          <div className="mx-auto max-w-6xl px-6 pt-28">
            <Reveal>
              <Eyebrow>{science.eyebrow}</Eyebrow>
            </Reveal>
            <Reveal delay={0.08}>
              <h2 className="mt-5 max-w-3xl text-3xl font-semibold tracking-tight sm:text-4xl">
                {science.title}
              </h2>
            </Reveal>
            <Reveal delay={0.16}>
              <p className="mt-5 max-w-2xl text-lg leading-relaxed text-muted">
                {science.sub}
              </p>
            </Reveal>
          </div>
          <HemeApproach />
        </section>

        {/* Values — how we work */}
        <section className="border-t border-border py-28">
          <div className="mx-auto max-w-4xl px-6">
            <Reveal>
              <Eyebrow>{values.eyebrow}</Eyebrow>
            </Reveal>
            <Reveal delay={0.08}>
              <blockquote className="mt-6 border-l-2 border-iron pl-6 text-2xl font-medium leading-snug tracking-tight text-fg sm:text-3xl">
                {values.title}
              </blockquote>
            </Reveal>
            <Reveal delay={0.16}>
              <p className="mt-8 max-w-2xl text-lg leading-relaxed text-muted">
                {values.body}
              </p>
            </Reveal>
          </div>
        </section>

        {/* Our first project — OpenAFR */}
        <section id="work" className="border-t border-border bg-surface/30 py-28">
          <div className="mx-auto max-w-6xl px-6">
            <Reveal>
              <Eyebrow>{projectSection.eyebrow}</Eyebrow>
            </Reveal>
            <Reveal delay={0.08}>
              <h2 className="mt-5 max-w-3xl text-3xl font-semibold tracking-tight sm:text-4xl">
                {projectSection.title}
              </h2>
            </Reveal>
            <Reveal delay={0.16}>
              <p className="mt-5 max-w-2xl text-lg leading-relaxed text-muted">
                {projectSection.intro}
              </p>
            </Reveal>

            <div className="mt-14 grid gap-5 sm:grid-cols-2 lg:grid-cols-4">
              {projectSection.stats.map((s) => (
                <Stat
                  key={s.label}
                  value={s.value}
                  label={s.label}
                  kind={s.kind}
                  suffix={"suffix" in s ? s.suffix : undefined}
                  tone={s.tone}
                />
              ))}
            </div>

            <Reveal delay={0.1}>
              <p className="mt-12 max-w-3xl text-pretty text-xl leading-relaxed text-fg">
                {projectSection.takeaway}
              </p>
            </Reveal>

            {/* Honest limits, kept with the project */}
            <div className="mt-14 rounded-2xl border border-border bg-bg/50 p-8">
              <h3 className="font-mono text-sm uppercase tracking-[0.15em] text-muted">
                {projectSection.limitsTitle}
              </h3>
              <ul className="mt-6 space-y-3">
                {projectSection.limits.map((pt) => (
                  <li key={pt} className="flex gap-4">
                    <span className="mt-1.5 h-1.5 w-1.5 shrink-0 rounded-full bg-faint" />
                    <span className="text-[15px] leading-relaxed text-muted">
                      {pt}
                    </span>
                  </li>
                ))}
              </ul>
            </div>

            <Reveal delay={0.1}>
              <a
                href={org.github}
                className="mt-10 inline-flex items-center gap-2 text-sm font-semibold text-iron transition-colors hover:text-fg"
              >
                Read the {project.name} code and every run write-up →
              </a>
            </Reveal>
          </div>
        </section>

        {/* Roadmap — where we're headed */}
        <section id="roadmap" className="border-t border-border py-28">
          <div className="mx-auto max-w-4xl px-6">
            <Reveal>
              <Eyebrow>{roadmap.eyebrow}</Eyebrow>
            </Reveal>
            <Reveal delay={0.08}>
              <h2 className="mt-5 text-3xl font-semibold tracking-tight sm:text-4xl">
                {roadmap.title}
              </h2>
            </Reveal>
            <ol className="mt-12 space-y-px">
              {roadmap.stages.map((st, i) => {
                const s = stateStyles[st.state];
                return (
                  <Reveal key={st.title} delay={i * 0.08} as="li">
                    <div className="relative flex gap-5 rounded-xl p-5 transition-colors hover:bg-surface/50">
                      <div className="flex flex-col items-center">
                        <span className={`mt-1.5 h-3 w-3 rounded-full ${s.dot}`} />
                        {i < roadmap.stages.length - 1 && (
                          <span className="mt-1 w-px flex-1 bg-border" />
                        )}
                      </div>
                      <div className="pb-6">
                        <div className="flex flex-wrap items-center gap-3">
                          <h3 className="text-lg font-semibold text-fg">
                            {st.title}
                          </h3>
                          <span
                            className={`font-mono text-[11px] uppercase tracking-wider ${s.text}`}
                          >
                            {s.label}
                          </span>
                        </div>
                        <p className="mt-2 max-w-xl text-[15px] leading-relaxed text-muted">
                          {st.body}
                        </p>
                      </div>
                    </div>
                  </Reveal>
                );
              })}
            </ol>
          </div>
        </section>

        {/* Contact / CTA */}
        <section
          id="contact"
          className="relative overflow-hidden border-t border-border py-32"
        >
          <div className="glow-iron pointer-events-none absolute inset-0" />
          <div className="relative mx-auto max-w-3xl px-6 text-center">
            <Reveal>
              <h2 className="text-balance text-3xl font-semibold tracking-tight sm:text-5xl">
                {cta.title}
              </h2>
            </Reveal>
            <Reveal delay={0.1}>
              <p className="mx-auto mt-6 max-w-2xl text-lg leading-relaxed text-muted">
                {cta.body}
              </p>
            </Reveal>
            <Reveal delay={0.2}>
              <div className="mt-10 flex flex-wrap justify-center gap-4">
                <a
                  href={`mailto:${org.email}`}
                  className="rounded-full bg-iron px-6 py-3 text-sm font-semibold text-bg transition-transform hover:-translate-y-0.5"
                >
                  Email us
                </a>
                <a
                  href={org.github}
                  className="rounded-full border border-border px-6 py-3 text-sm font-semibold text-fg transition-colors hover:border-iron hover:text-iron"
                >
                  Explore our work
                </a>
              </div>
            </Reveal>
          </div>
        </section>
      </main>

      <footer className="border-t border-border py-10">
        <div className="mx-auto flex max-w-6xl flex-col items-center justify-between gap-4 px-6 text-sm text-faint sm:flex-row">
          <span>
            {org.legal} · {org.mission}
          </span>
          <span className="font-mono text-xs">
            Open science · a rank is a hypothesis, never a hit
          </span>
        </div>
      </footer>
    </>
  );
}

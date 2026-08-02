/*
 * Single source of truth for site copy.
 * Edit here to change headings, stats, and CTAs — no component edits needed.
 *
 * The site leads with the FOUNDATION's mission (fungal research; a world without
 * antifungal resistance). OpenAFR — the screening pipeline — is presented as our
 * first project, not the headline.
 */

export const org = {
  name: "The Antifungal Resistance Project",
  short: "ARP", // compact mark for the nav
  legal: "The Antifungal Resistance Project",
  domain: "antifungalresistance.org", // TODO(you): register / confirm the domain
  mission: "For a world without antifungal resistance.",
  tagline: "A nonprofit for fungal research and open antifungal science.",
  email: "hello@antifungalresistance.org", // TODO(you): confirm the real inbox
  github: "https://github.com/", // TODO(you): repo / org URL
};

// The tool, referenced as a project throughout.
export const project = {
  name: "OpenAFR",
  blurb: "our open, self-validating antifungal screening pipeline",
};

export const hero = {
  eyebrow: org.tagline,
  title: "For a world without antifungal resistance.",
  sub: "Fungal disease is linked to an estimated 3.8 million deaths a year, the drugs that fight it are failing, and almost no one is working on it. We research the biology of fungal disease and build open, rigorously validated science to keep antifungal medicines working.",
};

export const problem = {
  eyebrow: "Why fungi, why now",
  title: "Fungal disease is rising, resistant, and overlooked.",
  points: [
    {
      k: "A hidden global killer",
      v: "Invasive fungal infections cause millions of deaths a year — often in people already fighting cancer, HIV, or ICU illness — yet they get a fraction of the attention bacteria do.",
    },
    {
      k: "Resistance is spreading",
      v: "Pathogens like the WHO critical-priority Candida auris now resist multiple antifungal classes at once and move through hospitals faster than we can respond.",
    },
    {
      k: "The pipeline is nearly empty",
      v: "Antifungal R&D is dramatically underfunded relative to antibacterial R&D. Few new drug classes are coming, and the ones we have are eroding.",
    },
  ],
};

export const approach = {
  eyebrow: "What we do",
  title: "Understand fungi. Build in the open. Hold every claim to proof.",
  pillars: [
    {
      k: "Research the biology",
      v: "Study how fungi cause disease and evade our drugs — starting with the enzymes and binding sites that antifungal resistance turns on.",
    },
    {
      k: "Build open tools",
      v: "Release reproducible, freely available software and data so any lab can build on the work, not just us.",
    },
    {
      k: "Bring people together",
      v: "Connect computational work with the wet labs, clinicians, and funders who can turn a hypothesis into a real treatment.",
    },
  ],
};

// Framing for the scroll-driven molecular animation (captions live in HemeApproach.tsx).
export const science = {
  eyebrow: "The science of resistance",
  title: "Resistance is a story about geometry.",
  sub: "Scroll to see how an antifungal drug actually works — and what really separates a molecule that kills fungus from one that only looks like it should.",
};

export const values = {
  eyebrow: "How we work",
  title: "We hold our own work to proof before we ask anyone to trust it.",
  body: "Computational biology produces hypotheses, not cures — and it is easy to publish a confident-looking result that means nothing. So we build the evidence in first: our tools must re-discover the drugs we already know work, on molecules they have never seen, before they are allowed to rank anything new. Every result, including the runs that fail, is written up in the open.",
};

// OpenAFR presented as the foundation's first project.
export const projectSection = {
  eyebrow: "Our first project",
  title: "OpenAFR — proving the method on Candida auris.",
  intro:
    "Our first open pipeline puts the approach to the test. On a held-out challenge frozen with a SHA-256 hash before any data was generated, it was asked to find 7 azole antifungals it had never seen, hidden among 348 look-alike decoys.",
  stats: [
    { label: "AUC — iron-approach distance", value: 0.794, kind: "decimal" as const, tone: "pass" as const },
    { label: "Enrichment @ 1%", value: 12.68, suffix: "×", kind: "decimal" as const, tone: "pass" as const },
    { label: "Permutation-test p", value: 0.0028, kind: "p" as const, tone: "pass" as const },
    { label: "AUC — docking affinity score", value: 0.471, kind: "decimal" as const, tone: "fail" as const },
  ],
  takeaway:
    "How closely a molecule brings a nitrogen to the heme iron discriminated real drugs from look-alikes far better than the docking program's own score — on the exact same poses. It's early, single-project evidence that the method works.",
  limitsTitle: "What it does not yet prove",
  limits: [
    "n = 7 held-out actives; the 95% confidence interval's lower bound still sits below the pass bar.",
    "Decoys are presumed inactive, not experimentally verified.",
    "One rigid receptor (C. albicans, not C. auris); no induced fit.",
    "This validates a ranking method, not a drug. A high rank is a hypothesis for a wet lab — never a hit.",
  ],
};

export const roadmap = {
  eyebrow: "Where we're headed",
  title: "From a proven method to a shared platform for antifungal discovery.",
  stages: [
    {
      state: "done" as const,
      title: "Prove the method",
      body: "Ship a self-validating pipeline whose ranking criterion passes a pre-registered, held-out test — with the failing runs published too.",
    },
    {
      state: "next" as const,
      title: "Broaden the science",
      body: "Extend beyond one enzyme and one pathogen, and turn ranked candidates into clear, reproducible reports a collaborator can act on.",
    },
    {
      state: "blocked" as const,
      title: "Partner with wet labs",
      body: "Put top candidates on real fungus. This is the critical dependency — nothing here is validated until someone does.",
    },
    {
      state: "next" as const,
      title: "Open the platform",
      body: "Grow into shared, open infrastructure and data that any antifungal researcher can use — so progress compounds across the whole field.",
    },
  ],
};

export const cta = {
  title: "Help us build a world without antifungal resistance.",
  body: "Whether you run a wet lab, study fungal biology, fund neglected disease work, or want to contribute to open science — there's a place for you here. Let's talk.",
};

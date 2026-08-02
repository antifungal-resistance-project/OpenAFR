# OpenAFR — foundation site

The public site for the OpenAFR foundation / PBC: what antifungal resistance is,
what we're building, and why the pipeline has to prove itself before it ranks
anything new.

Built to deploy on Vercel.

## Stack

- **Next.js 16** (App Router, Turbopack) + **React 19**
- **Tailwind CSS v4** (design tokens in `src/app/globals.css`)
- **Motion 12** (`motion/react`) for scroll-driven animation

## Develop

```bash
npm install
npm run dev      # http://localhost:3000
npm run build    # production build (fully static)
npm run lint
```

## Where to edit

- **All copy** lives in `src/content/site.ts` — headings, stats, roadmap, CTAs.
  Change wording there, not in components.
- **Colors / type / spacing tokens**: `src/app/globals.css` (`@theme`).
  `iron` = heme iron (primary accent), `nitrogen` = candidate molecule,
  `pass` = validated, `fail` = worse-than-random.
- **Sections & layout**: `src/app/page.tsx`.
- **The scroll animation**: `src/components/HemeApproach.tsx` — a candidate
  molecule's nitrogen approaching the heme iron as you scroll, with the N–Fe
  distance ticking to 2.63 Å. Geometry is derived from one scroll-progress value.
- Reusable pieces: `Reveal` (scroll fade-in), `Stat` (count-up), `Nav`.

Everything respects `prefers-reduced-motion`.

The site leads with the foundation's mission; **OpenAFR** appears as our first
project (the "Our work" section), not the headline.

## TODO before launch

- Pick the real foundation / PBC name and set `org.name`, `org.legal`,
  `org.domain` in `src/content/site.ts` (currently the placeholder "Mycora").
- Set `org.email` and `org.github`.
- Replace the default favicon in `src/app/favicon.ico`.
- Add real Open Graph / social images.

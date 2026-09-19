# Pre-rendered pocket images for the resistance weather report

The weather page (`openafr/weather.py`) references a wild-type-vs-mutant CYP51 pocket
image per alert at `pockets/<token>.png`. PyMOL is out of the offline verdict path
(see `openafr/structural.py`), so these PNGs are **pre-rendered out-of-band** and
committed here as repo assets; the page and the CI workflow only reference/copy them,
they never run PyMOL.

## How to render one (out-of-band, needs the pinned structural env + PyMOL)

For a buildable variant token (e.g. `F126L`), run the same chain `structural.py` emits
as its `render_hint`:

```bash
python scripts/mutate_receptor.py <token> -o work/receptor_auris_<token>.pdb \
  && python scripts/make_pose_view.py <azole> --receptor work/receptor_auris_<token>.pdb \
  && bash scripts/render_views.sh work/views/<token>
```

Then save the resulting image as `docs/pockets/<token>.png`, where `<token>` is the
mutation with non-alphanumerics stripped (matches `openafr.weather.pocket_asset`):
`F126L` -> `F126L.png`, `TR34/L98H` -> `TR34L98H.png`.

A missing PNG degrades gracefully: the page hides the broken image (`onerror`) and shows
the text verdict, so the feed is always valid even before any pocket has been rendered.

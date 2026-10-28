# OpenAFR — agent guide

## Project tracking

- **[GitHub Project "OpenAFR"](https://github.com/users/jajjer/projects/1)** is the operational
  board. Each active bet is an issue with a `Rank` + `Status` (Backlog / Next / In progress /
  Parked / Done) and a one-line `Why`.
- **[ROADMAP.md](ROADMAP.md)** is the narrative north-star doc: it tells the *why* behind each
  ranked bet and links to its issue. [VISION.md](VISION.md) = why the project exists;
  [TODOS.md](TODOS.md) = archive of deferred/finished work with cold-pickup context.
- **The loop:** pick the top `Next` / highest-rank live bet, ship **one PR**, then re-rank —
  move the item on the board and refresh ROADMAP's prose. One source of truth per fact: rank and
  status live on the board, the *why* lives in ROADMAP.
- **Every PR links its issue** (`Closes #<n>` in the body) so merging auto-moves the board item
  to Done. A CI job (`.github/workflows/pr-links-issue.yml`) nudges if you forget; opt out with
  the `no-issue` label for PRs that intentionally close nothing.
- **When you open a new roadmap issue, add it to the board** in the same breath:
  `gh project item-add 1 --owner jajjer --url <issue-url>` (board = user project #1, owner
  `jajjer`; set its `Rank`/`Status`/`Why` with `gh project item-edit`). This is how new work
  reaches the board — GitHub's UI auto-add workflow is not used.
- Not a sprint tracker. No scrum ceremony — this is a solo-plus-agents workflow; keep it light.

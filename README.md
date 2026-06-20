# OSCARLOCATOR Simulator — GitHub Pages build

A single-file build of the OSCARLOCATOR simulator for hosting on GitHub Pages.
Live AMSAT element data is fetched through **your** Cloudflare Worker, whose URL
is set in a config constant and is **never shown in the page interface**.

## What's here
- `index.html` — the deployable page (self-contained: engine, map data, and a
  bundled AMSAT snapshot are all inlined).
- `proxy/` — the Cloudflare Worker that re-serves the AMSAT GP bulletin with a
  CORS header (`worker.js`, `wrangler.toml`, `DEPLOY.md`).
- `simulator.ghpages.template.html` — the editable source (before inlining), if
  you want to rebuild.

## 1. Deploy the Worker
The AMSAT bulletin (`newark192.amsat.org`) is served without a CORS header, so a
browser cannot read it directly. The Worker in `proxy/` fixes that. See
`proxy/DEPLOY.md`; in short:

```
cd proxy
npx wrangler deploy
```

Note the deployed URL, e.g. `https://oscarlocator-pwa.YOURNAME.workers.dev`.

## 2. Set your Worker URL (one line)
Open `index.html`, find this block near the top of the script (search for
`CONFIG`) and replace the placeholder with your Worker URL:

```js
// ===================================================================
// CONFIG -- set this to your Cloudflare Worker URL before deploying.
// It is read only here and is never shown in the page UI.
const PROXY_URL = "https://YOUR-WORKER.workers.dev";
// ===================================================================
```

The URL lives only in this constant. It is not rendered anywhere in the page,
not placed in any input field, and not written to the DOM — so visitors don't
see it. (It is, of course, still visible to anyone who reads the page source or
watches network requests; a Worker URL is not a secret, but it is kept out of
the visible interface as requested. Lock the Worker down by origin in
`proxy/worker.js` if you want to restrict who can call it.)

## 3. Publish to GitHub Pages (custom subdomain)

This folder includes a `CNAME` file set to **oscarlocator.n8hm.radio**. To serve
the site from that subdomain:

1. Push this folder's contents (including `CNAME`) to a GitHub repo.
2. **Settings → Pages → Build and deployment → Source: Deploy from a branch**,
   pick your branch and `/ (root)`, save.
3. In your `n8hm.radio` DNS, add a CNAME record:
   `oscarlocator  CNAME  YOURUSERNAME.github.io`  (bare github.io host, not the repo URL).
4. Back in **Settings → Pages → Custom domain**, confirm `oscarlocator.n8hm.radio`
   is shown (the `CNAME` file sets this automatically), then tick **Enforce HTTPS**
   once DNS propagates. HTTPS is required for the geolocation auto-detect to work.

The site will be live at `https://oscarlocator.n8hm.radio/`.

The bundled Worker (`proxy/worker.js`) is **locked to this subdomain**: it only
returns the AMSAT bulletin to requests from `https://oscarlocator.n8hm.radio`,
`https://n8hm.radio`, and `http://localhost:8000`. Edit `ALLOW_ORIGINS` in
`proxy/worker.js` to change that list (or set it to `["*"]` to allow any origin),
then redeploy the Worker.

## 3b. Generic GitHub Pages (no custom domain)
1. Create a repo and add `index.html` (and `proxy/` if you like) to it.
2. Repo **Settings → Pages → Build and deployment → Source: Deploy from a
   branch**, pick your branch and `/ (root)`, save.
3. Your site appears at `https://YOURNAME.github.io/REPO/`.

## Behaviour
- On load, the page reads `PROXY_URL`; if it's been set (not the placeholder),
  it automatically fetches the current AMSAT catalog through your Worker and
  fills the satellite picker. **Refresh** re-fetches on demand.
- If the Worker is unreachable, the page falls back to the **bundled snapshot**
  baked into `index.html` so the picker still works, and shows a brief note.
- Everything else (projections, range circle, footprint, live tracking, the
  reference-orbit table, Maidenhead grid entry) works exactly as in the main
  build.

## Rebuilding index.html
If you edit `simulator.ghpages.template.html`, re-inline the assets (engine, land
GeoJSON, AMSAT snapshot) the same way the main build does, writing the result to
`index.html`. The three placeholders are the inline engine `<script>`, `__LAND__`,
and `__AMSAT__`.

## Refreshing the bundled snapshot
The fallback catalog is a point-in-time copy. To refresh it, fetch the bulletin
(through your Worker or any CORS proxy), slim it to the fields the sim uses
(`AMSAT_NAME, OBJECT_NAME, EPOCH, INCLINATION, RA_OF_ASC_NODE, ECCENTRICITY,
ARG_OF_PERICENTER, MEAN_ANOMALY, MEAN_MOTION, BSTAR`), and replace the
`__AMSAT__` block / `#amsat-data` script contents.

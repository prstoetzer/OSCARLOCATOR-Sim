# Deploying the AMSAT CORS proxy (Cloudflare Worker)

The AMSAT server doesn't send CORS headers, so a browser can't read its bulletin
directly. This tiny Worker fetches it server-side and re-serves it with the
right header. It's free (Cloudflare's free tier covers ~100k requests/day) and
stateless.

You only need to do this once. Two ways:

---

## Option A — Dashboard (no tools, ~3 minutes)

1. Sign up / log in at https://dash.cloudflare.com → **Workers & Pages**.
2. **Create application** → **Create Worker**. Give it a name, e.g.
   `oscarlocator-amsat-proxy`. Click **Deploy** (the placeholder code is fine
   for now).
3. Click **Edit code**. Delete everything in the editor and paste the entire
   contents of `worker.js` from this folder.
4. Click **Deploy**.
5. Copy your Worker URL — it looks like:
   `https://oscarlocator-amsat-proxy.YOURNAME.workers.dev`
6. Open it in a browser. You should see the AMSAT JSON. Done.

---

## Option B — Command line (Wrangler)

Requires Node.js.

```bash
cd proxy
npm install -g wrangler        # if you don't have it
wrangler login                 # opens a browser to authorize
wrangler deploy                # uses wrangler.toml + worker.js
```

Wrangler prints the deployed URL on success.

---

## Wire it into the app

1. Open `app.js` (in the parent folder) — or `src/app.jsx` if you recompile.
2. Find this line near the top of the config block:

   ```js
   const PROXY_URL = "";
   ```

3. Paste your Worker URL between the quotes:

   ```js
   const PROXY_URL = "https://oscarlocator-amsat-proxy.YOURNAME.workers.dev";
   ```

   If you edit `src/app.jsx`, recompile with:
   `npx babel src/app.jsx --presets @babel/preset-react -o app.js`
   (using classic runtime). If you edit `app.js` directly, just change the same
   line — it's plain JS in there too.

4. Re-host the app. "Fetch AMSAT live" will now use your proxy first, fall back
   to a public proxy, then to a direct attempt.

---

## Optional: lock the proxy to your site

In `worker.js`, change:

```js
const ALLOW_ORIGIN = "*";
```

to your app's exact origin, e.g.:

```js
const ALLOW_ORIGIN = "https://YOURNAME.github.io";
```

Then redeploy. This stops other sites from using your Worker. For read-only
public data it's not strictly necessary, but it's tidy.

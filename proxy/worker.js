/**
 * OSCARLOCATOR CORS proxy — Cloudflare Worker
 * --------------------------------------------------------------------------
 * Fetches the AMSAT daily GP-element bulletin server-side (where CORS does
 * not apply) and re-serves it with an Access-Control-Allow-Origin header so a
 * browser-based app can read it.
 *
 * The upstream URL is fixed, so this proxy can ONLY ever return the AMSAT
 * bulletin — it cannot be abused as an open proxy for arbitrary URLs.
 *
 * Deploy: see DEPLOY.md in this folder.
 */

const UPSTREAM = "https://newark192.amsat.org/gpdata/current/daily-bulletin.json";

// Lock the proxy to your own site(s). Requests from an origin in this list get
// that origin echoed back in Access-Control-Allow-Origin; anything else is
// refused at the CORS layer (the browser blocks the response). Add/remove
// origins as needed — each must be an exact scheme+host (no trailing slash/path).
// Set ALLOW_ORIGINS = ["*"] to allow any origin (open read-only access).
const ALLOW_ORIGINS = [
  "https://oscarlocator.n8hm.radio",
  "https://n8hm.radio",
  "http://localhost:8000",   // local testing; remove if you don't need it
];

// Edge cache lifetime in seconds. The AMSAT bulletin updates daily; caching for
// an hour keeps you well within that while being kind to their server.
const CACHE_SECONDS = 3600;

// Resolve the Access-Control-Allow-Origin value for this request:
//  - if ALLOW_ORIGINS contains "*", return "*"
//  - else if the request's Origin is allow-listed, echo it back
//  - else return null (no CORS header -> browser blocks cross-origin reads)
function allowedOrigin(request) {
  if (ALLOW_ORIGINS.includes("*")) return "*";
  const origin = request.headers.get("Origin");
  return origin && ALLOW_ORIGINS.includes(origin) ? origin : null;
}

function corsHeaders(allowOrigin, extra = {}) {
  const h = {
    "Access-Control-Allow-Methods": "GET, OPTIONS",
    "Access-Control-Allow-Headers": "Content-Type",
    "Access-Control-Max-Age": "86400",
    "Vary": "Origin",
    ...extra,
  };
  if (allowOrigin) h["Access-Control-Allow-Origin"] = allowOrigin;
  return h;
}

export default {
  async fetch(request) {
    const allowOrigin = allowedOrigin(request);

    // Preflight
    if (request.method === "OPTIONS") {
      return new Response(null, { status: 204, headers: corsHeaders(allowOrigin) });
    }
    if (request.method !== "GET") {
      return new Response("Method Not Allowed", {
        status: 405,
        headers: corsHeaders(allowOrigin, { "Content-Type": "text/plain" }),
      });
    }

    try {
      const upstream = await fetch(UPSTREAM, {
        // Cache at Cloudflare's edge so repeated app loads don't hammer AMSAT.
        cf: { cacheTtl: CACHE_SECONDS, cacheEverything: true },
        headers: { "Accept": "application/json" },
      });

      if (!upstream.ok) {
        return new Response(
          JSON.stringify({ error: "upstream", status: upstream.status }),
          { status: 502, headers: corsHeaders(allowOrigin, { "Content-Type": "application/json" }) }
        );
      }

      // Stream the body through, attach CORS + cache headers.
      const body = await upstream.text();
      return new Response(body, {
        status: 200,
        headers: corsHeaders(allowOrigin, {
          "Content-Type": "application/json; charset=utf-8",
          "Cache-Control": `public, max-age=${CACHE_SECONDS}`,
        }),
      });
    } catch (err) {
      return new Response(
        JSON.stringify({ error: "fetch_failed", message: String(err) }),
        { status: 502, headers: corsHeaders(allowOrigin, { "Content-Type": "application/json" }) }
      );
    }
  },
};

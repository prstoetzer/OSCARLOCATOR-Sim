(function(){
'use strict';
/* =========================================================================
 * sgp4.js — pure-JavaScript SGP4/SDP4 propagator (WGS72)
 * -------------------------------------------------------------------------
 * A faithful port of the Vallado "Revisiting Spacetrack Report #3"
 * (AIAA-2006-6753) reference implementation, including the deep-space (SDP4)
 * terms. No dependencies. Mirrors the math used by OrbitDeck's sgp4_lite.py
 * (both derive from the same Vallado reference), so positions agree with the
 * desktop app and with the reference `sgp4` package to well under a km for LEO.
 *
 * Public API:
 *   const sat = sgp4Init({ inclo, nodeo, ecco, argpo, mo, no_kozai, bstar,
 *                          epochjd });            // angles in RADIANS, no_kozai rev->rad/min
 *   const { r, v, error } = sgp4(sat, tsince);    // tsince in MINUTES from epoch
 *     r,v are km and km/s in TEME.
 *
 * Convenience:
 *   sgp4FromOMM(omm)   -> sat   (accepts AMSAT/CelesTrak OMM-JSON field names)
 *   propagateAt(sat, jsDate) -> { r, v, error }
 * ========================================================================= */

const pi = Math.PI;
const twopi = 2.0 * pi;
const deg2rad = pi / 180.0;
const x2o3 = 2.0 / 3.0;

// WGS72 gravity constants (the SGP4 default)
const mu = 398600.8;            // km^3/s^2
const radiusearthkm = 6378.135; // km
const xke = 60.0 / Math.sqrt((radiusearthkm * radiusearthkm * radiusearthkm) / mu);
const tumin = 1.0 / xke;
const j2 = 0.001082616;
const j3 = -0.00000253881;
const j4 = -0.00000165597;
const j3oj2 = j3 / j2;

function gstime(jdut1) {
  let tut1 = (jdut1 - 2451545.0) / 36525.0;
  let temp =
    -6.2e-6 * tut1 * tut1 * tut1 +
    0.093104 * tut1 * tut1 +
    (876600.0 * 3600 + 8640184.812866) * tut1 +
    67310.54841;
  temp = ((temp * deg2rad) / 240.0) % twopi;
  if (temp < 0.0) temp += twopi;
  return temp;
}

/* ---- deep-space initialization (dscom, dpper, dsinit, dspace) ---------- */
// These follow Vallado's structure closely. Kept compact but complete.

function dpper(e3, ee2, peo, pgho, pho, pinco, plo, se2, se3, sgh2, sgh3, sgh4,
  sh2, sh3, si2, si3, sl2, sl3, sl4, t, xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3,
  xl2, xl3, xl4, zmol, zmos, init, op) {
  const zns = 1.19459e-5, zes = 0.01675, znl = 1.5835218e-4, zel = 0.05490;
  let zm, zf, sinzf, f2, f3, ses, sis, sls, sghs, shll, sel, sil, sll, sghl, shl;
  let pe = op.ep, pinc = op.inclp, pl = op.mp, pgh = op.argpp, ph = op.nodep;

  zm = zmos + zns * t;
  if (init === "y") zm = zmos;
  zf = zm + 2.0 * zes * Math.sin(zm);
  sinzf = Math.sin(zf);
  f2 = 0.5 * sinzf * sinzf - 0.25;
  f3 = -0.5 * sinzf * Math.cos(zf);
  ses = se2 * f2 + se3 * f3;
  sis = si2 * f2 + si3 * f3;
  sls = sl2 * f2 + sl3 * f3 + sl4 * sinzf;
  sghs = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
  shll = sh2 * f2 + sh3 * f3;

  zm = zmol + znl * t;
  if (init === "y") zm = zmol;
  zf = zm + 2.0 * zel * Math.sin(zm);
  sinzf = Math.sin(zf);
  f2 = 0.5 * sinzf * sinzf - 0.25;
  f3 = -0.5 * sinzf * Math.cos(zf);
  sel = ee2 * f2 + e3 * f3;
  sil = xi2 * f2 + xi3 * f3;
  sll = xl2 * f2 + xl3 * f3 + xl4 * sinzf;
  sghl = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
  shl = xh2 * f2 + xh3 * f3;

  pe = ses + sel; pinc = sis + sil; pl = sls + sll;
  pgh = sghs + sghl; ph = shll + shl;

  if (init === "n") {
    pe = pe - peo; pinc = pinc - pinco; pl = pl - plo;
    pgh = pgh - pgho; ph = ph - pho;
    op.inclp += pinc; op.ep += pe;
    const sinip = Math.sin(op.inclp), cosip = Math.cos(op.inclp);

    if (op.inclp >= 0.2) {
      ph = ph / sinip; pgh = pgh - cosip * ph;
      op.argpp += pgh; op.nodep += ph; op.mp += pl;
    } else {
      const sinop = Math.sin(op.nodep), cosop = Math.cos(op.nodep);
      let alfdp = sinip * sinop, betdp = sinip * cosop;
      const dalf = ph * cosop + pinc * cosip * sinop;
      const dbet = -ph * sinop + pinc * cosip * cosop;
      alfdp += dalf; betdp += dbet;
      op.nodep = op.nodep % twopi;
      if (op.nodep < 0.0) op.nodep += twopi;
      let xls = op.mp + op.argpp + cosip * op.nodep;
      let dls = pl + pgh - pinc * op.nodep * sinip;
      xls += dls;
      const xnoh = op.nodep;
      op.nodep = Math.atan2(alfdp, betdp);
      if (op.nodep < 0.0) op.nodep += twopi;
      if (Math.abs(xnoh - op.nodep) > pi) {
        if (op.nodep < xnoh) op.nodep += twopi; else op.nodep -= twopi;
      }
      op.mp += pl;
      op.argpp = xls - op.mp - cosip * op.nodep;
    }
  }
}

function dscom(epoch, ep, argpp, tc, inclp, nodep, np) {
  const zes = 0.01675, zel = 0.05490, c1ss = 2.9864797e-6, c1l = 4.7968065e-7;
  const zsinis = 0.39785416, zcosis = 0.91744867, zcosgs = 0.1945905, zsings = -0.98088458;
  const o = {};
  o.nm = np; o.em = ep; o.snodm = Math.sin(nodep); o.cnodm = Math.cos(nodep);
  o.sinomm = Math.sin(argpp); o.cosomm = Math.cos(argpp);
  o.sinim = Math.sin(inclp); o.cosim = Math.cos(inclp);
  o.emsq = o.em * o.em; const betasq = 1.0 - o.emsq; o.rtemsq = Math.sqrt(betasq);
  o.peo = 0; o.pinco = 0; o.plo = 0; o.pgho = 0; o.pho = 0;
  o.day = epoch + 18261.5 + tc / 1440.0;
  const xnodce = (4.5236020 - 9.2422029e-4 * o.day) % twopi;
  const stem = Math.sin(xnodce), ctem = Math.cos(xnodce);
  const zcosil = 0.91375164 - 0.03568096 * ctem;
  const zsinil = Math.sqrt(1.0 - zcosil * zcosil);
  const zsinhl = (0.089683511 * stem) / zsinil;
  const zcoshl = Math.sqrt(1.0 - zsinhl * zsinhl);
  o.gam = 5.8351514 + 0.0019443680 * o.day;
  let zx = 0.39785416 * stem / zsinil;
  const zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
  zx = Math.atan2(zx, zy); zx = o.gam + zx - xnodce;
  const zcosgl = Math.cos(zx), zsingl = Math.sin(zx);

  let zcosg = zcosgs, zsing = zsings, zcosi = zcosis, zsini = zsinis;
  let zcosh = o.cnodm, zsinh = o.snodm, cc = c1ss, xnoi = 1.0 / o.nm;

  o.e3 = 0; o.ee2 = 0; o.se2 = 0; o.se3 = 0; o.sgh2 = 0; o.sgh3 = 0; o.sgh4 = 0;
  o.sh2 = 0; o.sh3 = 0; o.si2 = 0; o.si3 = 0; o.sl2 = 0; o.sl3 = 0; o.sl4 = 0;
  o.xgh2 = 0; o.xgh3 = 0; o.xgh4 = 0; o.xh2 = 0; o.xh3 = 0; o.xi2 = 0; o.xi3 = 0;
  o.xl2 = 0; o.xl3 = 0; o.xl4 = 0; o.zmol = 0; o.zmos = 0;

  for (let lsflg = 1; lsflg <= 2; lsflg++) {
    const a1 = zcosg * zcosh + zsing * zcosi * zsinh;
    const a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
    const a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
    const a8 = zsing * zsini;
    const a9 = zsing * zsinh + zcosg * zcosi * zcosh;
    const a10 = zcosg * zsini;
    const a2 = o.cosim * a7 + o.sinim * a8;
    const a4 = o.cosim * a9 + o.sinim * a10;
    const a5 = -o.sinim * a7 + o.cosim * a8;
    const a6 = -o.sinim * a9 + o.cosim * a10;
    const x1 = a1 * o.cosomm + a2 * o.sinomm;
    const x2 = a3 * o.cosomm + a4 * o.sinomm;
    const x3 = -a1 * o.sinomm + a2 * o.cosomm;
    const x4 = -a3 * o.sinomm + a4 * o.cosomm;
    const x5 = a5 * o.sinomm, x6 = a6 * o.sinomm;
    const x7 = a5 * o.cosomm, x8 = a6 * o.cosomm;
    o.z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
    o.z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
    o.z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
    o.z1 = 3.0 * (a1 * a1 + a2 * a2) + o.z31 * o.emsq;
    o.z2 = 6.0 * (a1 * a3 + a2 * a4) + o.z32 * o.emsq;
    o.z3 = 3.0 * (a3 * a3 + a4 * a4) + o.z33 * o.emsq;
    o.z11 = -6.0 * a1 * a5 + o.emsq * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
    o.z12 = -6.0 * (a1 * a6 + a3 * a5) + o.emsq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
    o.z13 = -6.0 * a3 * a6 + o.emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
    o.z21 = 6.0 * a2 * a5 + o.emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
    o.z22 = 6.0 * (a4 * a5 + a2 * a6) + o.emsq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
    o.z23 = 6.0 * a4 * a6 + o.emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
    o.z1 = o.z1 + o.z1 + betasq * o.z31;
    o.z2 = o.z2 + o.z2 + betasq * o.z32;
    o.z3 = o.z3 + o.z3 + betasq * o.z33;
    o.s3 = cc * xnoi;
    o.s2 = (-0.5 * o.s3) / o.rtemsq;
    o.s4 = o.s3 * o.rtemsq;
    o.s1 = -15.0 * o.em * o.s4;
    o.s5 = x1 * x3 + x2 * x4;
    o.s6 = x2 * x3 + x1 * x4;
    o.s7 = x2 * x4 - x1 * x3;
    if (lsflg === 1) {
      o.ss1 = o.s1; o.ss2 = o.s2; o.ss3 = o.s3; o.ss4 = o.s4; o.ss5 = o.s5;
      o.ss6 = o.s6; o.ss7 = o.s7; o.sz1 = o.z1; o.sz2 = o.z2; o.sz3 = o.z3;
      o.sz11 = o.z11; o.sz12 = o.z12; o.sz13 = o.z13; o.sz21 = o.z21; o.sz22 = o.z22;
      o.sz23 = o.z23; o.sz31 = o.z31; o.sz32 = o.z32; o.sz33 = o.z33;
      zcosg = zcosgl; zsing = zsingl; zcosi = zcosil; zsini = zsinil;
      zcosh = zcoshl * o.cnodm + zsinhl * o.snodm;
      zsinh = o.snodm * zcoshl - o.cnodm * zsinhl;
      cc = c1l;
    }
  }
  o.zmol = (4.7199672 + 0.22997150 * o.day - o.gam) % twopi;
  o.zmos = (6.2565837 + 0.017201977 * o.day) % twopi;
  if (o.zmol < 0) o.zmol += twopi;
  if (o.zmos < 0) o.zmos += twopi;

  o.se2 = 2.0 * o.ss1 * o.ss6; o.se3 = 2.0 * o.ss1 * o.ss7;
  o.si2 = 2.0 * o.ss2 * o.sz12; o.si3 = 2.0 * o.ss2 * (o.sz13 - o.sz11);
  o.sl2 = -2.0 * o.ss3 * o.sz2; o.sl3 = -2.0 * o.ss3 * (o.sz3 - o.sz1);
  o.sl4 = -2.0 * o.ss3 * (-21.0 - 9.0 * o.emsq) * zes;
  o.sgh2 = 2.0 * o.ss4 * o.sz32; o.sgh3 = 2.0 * o.ss4 * (o.sz33 - o.sz31);
  o.sgh4 = -18.0 * o.ss4 * zes; o.sh2 = -2.0 * o.ss2 * o.sz22;
  o.sh3 = -2.0 * o.ss2 * (o.sz23 - o.sz21);
  o.ee2 = 2.0 * o.s1 * o.s6; o.e3 = 2.0 * o.s1 * o.s7;
  o.xi2 = 2.0 * o.s2 * o.z12; o.xi3 = 2.0 * o.s2 * (o.z13 - o.z11);
  o.xl2 = -2.0 * o.s3 * o.z2; o.xl3 = -2.0 * o.s3 * (o.z3 - o.z1);
  o.xl4 = -2.0 * o.s3 * (-21.0 - 9.0 * o.emsq) * zel;
  o.xgh2 = 2.0 * o.s4 * o.z32; o.xgh3 = 2.0 * o.s4 * (o.z33 - o.z31);
  o.xgh4 = -18.0 * o.s4 * zel; o.xh2 = -2.0 * o.s2 * o.z22;
  o.xh3 = -2.0 * o.s2 * (o.z23 - o.z21);
  return o;
}

// dsinit: writes all secular-rate and resonance terms onto `d` (the deep
// struct that dspace later reads), using scalars passed in `args`.
function dsinit(d, tc, xpidot, args) {
  const q22 = 1.7891679e-6, q31 = 2.1460748e-6, q33 = 2.2123015e-7;
  const root22 = 1.7891679e-6, root44 = 7.3636953e-9, root54 = 2.1765803e-9;
  const rptim = 4.37526908801129966e-3, root32 = 3.7393792e-7, root52 = 1.1428639e-7;
  const znl = 1.5835218e-4, zns = 1.19459e-5;
  const s = {
    nm: args.no_unkozai, em: args.ecco, emsq: args.eccsq, inclm: args.inclo,
    argpm: 0.0, nodem: 0.0, mm: 0.0, t: 0.0,
  };
  let aonv;
  d.irez = 0;
  if (s.nm < 0.0052359877 && s.nm > 0.0034906585) d.irez = 1;
  if (s.nm >= 8.26e-3 && s.nm <= 9.24e-3 && s.em >= 0.5) d.irez = 2;

  // solar / lunar terms
  let ses = d.ss1 * zns * d.ss5;
  let sis = d.ss2 * zns * (d.sz11 + d.sz13);
  let sls = -zns * d.ss3 * (d.sz1 + d.sz3 - 14.0 - 6.0 * d.emsq);
  let sghs = d.ss4 * zns * (d.sz31 + d.sz33 - 6.0);
  let shs = -zns * d.ss2 * (d.sz21 + d.sz23);
  if (s.inclm < 5.2359877e-2 || s.inclm > pi - 5.2359877e-2) shs = 0.0;
  if (d.sinim !== 0.0) shs = shs / d.sinim;
  let sgs = sghs - d.cosim * shs;
  d.dedt = ses + d.s1 * znl * d.s5;
  d.didt = sis + d.s2 * znl * (d.z11 + d.z13);
  d.dmdt = sls - znl * d.s3 * (d.z1 + d.z3 - 14.0 - 6.0 * d.emsq);
  let sghl = d.s4 * znl * (d.z31 + d.z33 - 6.0);
  let shll = -znl * d.s2 * (d.z21 + d.z23);
  if (s.inclm < 5.2359877e-2 || s.inclm > pi - 5.2359877e-2) shll = 0.0;
  d.domdt = sgs + sghl;
  d.dnodt = shs;
  if (d.sinim !== 0.0) {
    d.domdt = d.domdt - d.cosim / d.sinim * shll;
    d.dnodt = d.dnodt + shll / d.sinim;
  }
  let dndt = 0.0;
  const theta = (args.gsto + tc * rptim) % twopi;
  // (init: args.t == 0, so the secular advances below are no-ops at epoch)
  s.em += d.dedt * args.t; s.inclm += d.didt * args.t;
  s.argpm += d.domdt * args.t; s.nodem += d.dnodt * args.t; s.mm += d.dmdt * args.t;

  // init resonance terms to zero so dspace sees defined values
  d.del1 = 0.0; d.del2 = 0.0; d.del3 = 0.0;
  d.d2201 = 0.0; d.d2211 = 0.0; d.d3210 = 0.0; d.d3222 = 0.0;
  d.d4410 = 0.0; d.d4422 = 0.0; d.d5220 = 0.0; d.d5232 = 0.0;
  d.d5421 = 0.0; d.d5433 = 0.0;
  d.xlamo = 0.0; d.xfact = 0.0; d.xli = 0.0; d.xni = 0.0; d.atime = 0.0;

  if (d.irez !== 0) {
    aonv = Math.pow(s.nm / xke, x2o3);
    if (d.irez === 2) {
      const cosisq = d.cosim * d.cosim;
      let emo = s.em, emsqo = s.emsq;
      s.em = args.ecco; s.emsq = args.eccsq;
      const eoc = s.em * s.emsq;
      let g201 = -0.306 - (s.em - 0.64) * 0.440;
      let g211, g310, g322, g410, g422, g520, g521, g532, g533;
      if (s.em <= 0.65) {
        g211 = 3.616 - 13.2470 * s.em + 16.2900 * s.emsq;
        g310 = -19.302 + 117.3900 * s.em - 228.4190 * s.emsq + 156.5910 * eoc;
        g322 = -18.9068 + 109.7927 * s.em - 214.6334 * s.emsq + 146.5816 * eoc;
        g410 = -41.122 + 242.6940 * s.em - 471.0940 * s.emsq + 313.9530 * eoc;
        g422 = -146.407 + 841.8800 * s.em - 1629.014 * s.emsq + 1083.435 * eoc;
        g520 = -532.114 + 3017.977 * s.em - 5740.032 * s.emsq + 3708.276 * eoc;
      } else {
        g211 = -72.099 + 331.819 * s.em - 508.738 * s.emsq + 266.724 * eoc;
        g310 = -346.844 + 1582.851 * s.em - 2415.925 * s.emsq + 1246.113 * eoc;
        g322 = -342.585 + 1554.908 * s.em - 2366.899 * s.emsq + 1215.972 * eoc;
        g410 = -1052.797 + 4758.686 * s.em - 7193.992 * s.emsq + 3651.957 * eoc;
        g422 = -3581.690 + 16178.110 * s.em - 24462.770 * s.emsq + 12422.520 * eoc;
        if (s.em > 0.715)
          g520 = -5149.66 + 29936.92 * s.em - 54087.36 * s.emsq + 31324.56 * eoc;
        else
          g520 = 1464.74 - 4664.75 * s.em + 3763.64 * s.emsq;
      }
      if (s.em < 0.7) {
        g533 = -919.22770 + 4988.6100 * s.em - 9064.7700 * s.emsq + 5542.21 * eoc;
        g521 = -822.71072 + 4568.6173 * s.em - 8491.4146 * s.emsq + 5337.524 * eoc;
        g532 = -853.66600 + 4690.2500 * s.em - 8624.7700 * s.emsq + 5341.4 * eoc;
      } else {
        g533 = -37995.780 + 161616.52 * s.em - 229838.20 * s.emsq + 109377.94 * eoc;
        g521 = -51752.104 + 218913.95 * s.em - 309468.16 * s.emsq + 146349.42 * eoc;
        g532 = -40023.880 + 170470.89 * s.em - 242699.48 * s.emsq + 115605.82 * eoc;
      }
      const sini2 = d.sinim * d.sinim;
      const f220 = 0.75 * (1.0 + 2.0 * d.cosim + cosisq);
      const f221 = 1.5 * sini2;
      const f321 = 1.875 * d.sinim * (1.0 - 2.0 * d.cosim - 3.0 * cosisq);
      const f322 = -1.875 * d.sinim * (1.0 + 2.0 * d.cosim - 3.0 * cosisq);
      const f441 = 35.0 * sini2 * f220;
      const f442 = 39.3750 * sini2 * sini2;
      const f522 = 9.84375 * d.sinim * (sini2 * (1.0 - 2.0 * d.cosim - 5.0 * cosisq) +
        0.33333333 * (-2.0 + 4.0 * d.cosim + 6.0 * cosisq));
      const f523 = d.sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * d.cosim + 10.0 * cosisq) +
        6.56250012 * (1.0 + 2.0 * d.cosim - 3.0 * cosisq));
      const f542 = 29.53125 * d.sinim * (2.0 - 8.0 * d.cosim + cosisq * (-12.0 + 8.0 * d.cosim + 10.0 * cosisq));
      const f543 = 29.53125 * d.sinim * (-2.0 - 8.0 * d.cosim + cosisq * (12.0 + 8.0 * d.cosim - 10.0 * cosisq));
      const xno2 = s.nm * s.nm;
      const ainv2 = aonv * aonv;
      let temp1 = 3.0 * xno2 * ainv2;
      let temp = temp1 * root22;
      d.d2201 = temp * f220 * g201;
      d.d2211 = temp * f221 * g211;
      temp1 = temp1 * aonv;
      temp = temp1 * root32;
      d.d3210 = temp * f321 * g310;
      d.d3222 = temp * f322 * g322;
      temp1 = temp1 * aonv;
      temp = 2.0 * temp1 * root44;
      d.d4410 = temp * f441 * g410;
      d.d4422 = temp * f442 * g422;
      temp1 = temp1 * aonv;
      temp = temp1 * root52;
      d.d5220 = temp * f522 * g520;
      d.d5232 = temp * f523 * g532;
      temp = 2.0 * temp1 * root54;
      d.d5421 = temp * f542 * g521;
      d.d5433 = temp * f543 * g533;
      d.xlamo = (args.mo + args.nodeo + args.nodeo - theta - theta) % twopi;
      d.xfact = args.mdot + d.dmdt + 2.0 * (args.nodedot + d.dnodt - rptim) - args.no_unkozai;
      s.em = emo; s.emsq = emsqo;
    }
    if (d.irez === 1) {
      const g200 = 1.0 + s.emsq * (-2.5 + 0.8125 * s.emsq);
      const g310 = 1.0 + 2.0 * s.emsq;
      const g300 = 1.0 + s.emsq * (-6.0 + 6.60937 * s.emsq);
      const f220 = 0.75 * (1.0 + d.cosim) * (1.0 + d.cosim);
      const f311 = 0.9375 * d.sinim * d.sinim * (1.0 + 3.0 * d.cosim) - 0.75 * (1.0 + d.cosim);
      let f330 = 1.0 + d.cosim;
      f330 = 1.875 * f330 * f330 * f330;
      d.del1 = 3.0 * s.nm * s.nm * aonv * aonv;
      d.del2 = 2.0 * d.del1 * f220 * g200 * q22;
      d.del3 = 3.0 * d.del1 * f330 * g300 * q33 * aonv;
      d.del1 = d.del1 * f311 * g310 * q31 * aonv;
      d.xlamo = (args.mo + args.nodeo + args.argpo - theta) % twopi;
      d.xfact = args.mdot + xpidot - rptim + d.dmdt + d.domdt + d.dnodt - args.no_unkozai;
    }
    d.xli = d.xlamo;
    d.xni = args.no_unkozai;
    d.atime = 0.0;
    d.nmInit = args.no_unkozai + dndt;
  }
}

function dspace(d, tc, op) {
  const fasx2 = 0.13130908, fasx4 = 2.8843198, fasx6 = 0.37448087;
  const g22 = 5.7686396, g32 = 0.95240898, g44 = 1.8014998, g52 = 1.0508330, g54 = 4.4108898;
  const rptim = 4.37526908801129966e-3, stepp = 720.0, stepn = -720.0, step2 = 259200.0;
  const s = op;
  let xndt = 0, xnddt = 0, xldot = 0, ft = 0;
  s.dndt = 0.0;
  const theta = (op.gsto + tc * rptim) % twopi;
  s.em += d.dedt * op.t; s.inclm += d.didt * op.t;
  s.argpm += d.domdt * op.t; s.nodem += d.dnodt * op.t; s.mm += d.dmdt * op.t;

  if (d.irez !== 0) {
    // integrator state persists on d across calls (Vallado satrec.atime/xli/xni)
    if (d.atime === 0.0 || op.t * d.atime <= 0.0 || Math.abs(op.t) < Math.abs(d.atime)) {
      d.atime = 0.0; d.xni = op.no_unkozai; d.xli = d.xlamo;
    }
    let delt;
    if (op.t > 0.0) delt = stepp; else delt = stepn;
    let iretn = 381;
    while (iretn === 381) {
      if (d.irez !== 2) {
        xndt = d.del1 * Math.sin(d.xli - fasx2) + d.del2 * Math.sin(2.0 * (d.xli - fasx4)) +
          d.del3 * Math.sin(3.0 * (d.xli - fasx6));
        xldot = d.xni + d.xfact;
        xnddt = d.del1 * Math.cos(d.xli - fasx2) + 2.0 * d.del2 * Math.cos(2.0 * (d.xli - fasx4)) +
          3.0 * d.del3 * Math.cos(3.0 * (d.xli - fasx6));
        xnddt = xnddt * xldot;
      } else {
        const xomi = op.argpo + op.argpdot * d.atime;
        const x2omi = xomi + xomi;
        const x2li = d.xli + d.xli;
        xndt = d.d2201 * Math.sin(x2omi + d.xli - g22) + d.d2211 * Math.sin(d.xli - g22) +
          d.d3210 * Math.sin(xomi + d.xli - g32) + d.d3222 * Math.sin(-xomi + d.xli - g32) +
          d.d4410 * Math.sin(x2omi + x2li - g44) + d.d4422 * Math.sin(x2li - g44) +
          d.d5220 * Math.sin(xomi + d.xli - g52) + d.d5232 * Math.sin(-xomi + d.xli - g52) +
          d.d5421 * Math.sin(xomi + x2li - g54) + d.d5433 * Math.sin(-xomi + x2li - g54);
        xldot = d.xni + d.xfact;
        xnddt = d.d2201 * Math.cos(x2omi + d.xli - g22) + d.d2211 * Math.cos(d.xli - g22) +
          d.d3210 * Math.cos(xomi + d.xli - g32) + d.d3222 * Math.cos(-xomi + d.xli - g32) +
          d.d5220 * Math.cos(xomi + d.xli - g52) + d.d5232 * Math.cos(-xomi + d.xli - g52) +
          2.0 * (d.d4410 * Math.cos(x2omi + x2li - g44) + d.d4422 * Math.cos(x2li - g44) +
            d.d5421 * Math.cos(xomi + x2li - g54) + d.d5433 * Math.cos(-xomi + x2li - g54));
        xnddt = xnddt * xldot;
      }
      if (Math.abs(op.t - d.atime) >= stepp) {
        iretn = 381;
      } else {
        ft = op.t - d.atime;
        iretn = 0;
      }
      if (iretn === 381) {
        d.xli = d.xli + xldot * delt + xndt * step2;
        d.xni = d.xni + xndt * delt + xnddt * step2;
        d.atime = d.atime + delt;
      }
    }
    s.nm = d.xni + xndt * ft + xnddt * ft * ft * 0.5;
    const xl = d.xli + xldot * ft + xndt * ft * ft * 0.5;
    if (d.irez !== 1) {
      s.mm = xl - 2.0 * s.nodem + 2.0 * theta;
      s.dndt = s.nm - op.no_unkozai;
    } else {
      s.mm = xl - s.nodem - s.argpm + theta;
      s.dndt = s.nm - op.no_unkozai;
    }
    s.nm = op.no_unkozai + s.dndt;
  }
}

/* ---- initialization ---------------------------------------------------- */
function sgp4Init(elem) {
  const s = {};
  s.error = 0;
  s.inclo = elem.inclo; s.nodeo = elem.nodeo; s.ecco = elem.ecco;
  s.argpo = elem.argpo; s.mo = elem.mo; s.bstar = elem.bstar || 0.0;
  s.no_kozai = elem.no_kozai; s.epochjd = elem.epochjd;

  s.t = 0.0;
  s.method = "n";
  const ss = 78.0 / radiusearthkm + 1.0;
  const qzms2t = Math.pow((120.0 - 78.0) / radiusearthkm, 4);

  // recover original mean motion (un-kozai)
  const eccsq = s.ecco * s.ecco;
  const omeosq = 1.0 - eccsq;
  const rteosq = Math.sqrt(omeosq);
  const cosio = Math.cos(s.inclo);
  const cosio2 = cosio * cosio;
  const ak = Math.pow(xke / s.no_kozai, x2o3);
  const d1 = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
  let del_ = d1 / (ak * ak);
  const adel = ak * (1.0 - del_ * del_ - del_ * (1.0 / 3.0 + 134.0 * del_ * del_ / 81.0));
  del_ = d1 / (adel * adel);
  const no = s.no_kozai / (1.0 + del_);
  s.no_unkozai = no;
  s.eccsq = eccsq;

  const ao = Math.pow(xke / no, x2o3);
  const sinio = Math.sin(s.inclo);
  const po = ao * omeosq;
  const con42 = 1.0 - 5.0 * cosio2;
  const con41 = -con42 - cosio2 - cosio2;
  s.con41 = con41;
  const ainv = 1.0 / ao;
  const posq = po * po;
  const rp = ao * (1.0 - s.ecco);

  s.gsto = gstime(s.epochjd);

  let sfour = ss;
  let qzms24 = qzms2t;
  const perige = (rp - 1.0) * radiusearthkm;
  if (perige < 156.0) {
    sfour = perige - 78.0;
    if (perige < 98.0) sfour = 20.0;
    qzms24 = Math.pow((120.0 - sfour) / radiusearthkm, 4.0);
    sfour = sfour / radiusearthkm + 1.0;
  }
  const pinvsq = 1.0 / posq;
  const tsi = 1.0 / (ao - sfour);
  s.eta = ao * s.ecco * tsi;
  const etasq = s.eta * s.eta;
  const eeta = s.ecco * s.eta;
  const psisq = Math.abs(1.0 - etasq);
  const coef = qzms24 * Math.pow(tsi, 4.0);
  const coef1 = coef / Math.pow(psisq, 3.5);
  const cc2 = coef1 * no * (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
    0.375 * j2 * tsi / psisq * con41 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
  s.cc1 = s.bstar * cc2;
  let cc3 = 0.0;
  if (s.ecco > 1.0e-4) cc3 = -2.0 * coef * tsi * j3oj2 * no * sinio / s.ecco;
  s.x1mth2 = 1.0 - cosio2;
  s.cc4 = 2.0 * no * coef1 * ao * omeosq *
    (s.eta * (2.0 + 0.5 * etasq) + s.ecco * (0.5 + 2.0 * etasq) -
      j2 * tsi / (ao * psisq) * (-3.0 * con41 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) +
        0.75 * s.x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) * Math.cos(2.0 * s.argpo)));
  s.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
  const cosio4 = cosio2 * cosio2;
  const temp1 = 1.5 * j2 * pinvsq * no;
  const temp2 = 0.5 * temp1 * j2 * pinvsq;
  const temp3 = -0.46875 * j4 * pinvsq * pinvsq * no;
  s.mdot = no + 0.5 * temp1 * rteosq * con41 + 0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
  s.argpdot = -0.5 * temp1 * con42 + 0.0625 * temp2 * (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
    temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
  const xhdot1 = -temp1 * cosio;
  s.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) + 2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
  s.xpidot = s.argpdot + s.nodedot;
  s.omgcof = s.bstar * cc3 * Math.cos(s.argpo);
  s.xmcof = 0.0;
  if (s.ecco > 1.0e-4) s.xmcof = -x2o3 * coef * s.bstar / eeta;
  s.nodecf = 3.5 * omeosq * xhdot1 * s.cc1;
  s.t2cof = 1.5 * s.cc1;
  if (Math.abs(cosio + 1.0) > 1.5e-12)
    s.xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
  else
    s.xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / 1.5e-12;
  s.aycof = -0.5 * j3oj2 * sinio;
  const delmotemp = 1.0 + s.eta * Math.cos(s.mo);
  s.delmo = delmotemp * delmotemp * delmotemp;
  s.sinmao = Math.sin(s.mo);
  s.x7thm1 = 7.0 * cosio2 - 1.0;

  // simplified-drag flag and higher-order secular drag coefficients
  s.isimp = 0;
  if (rp < (220.0 / radiusearthkm + 1.0)) s.isimp = 1;
  s.d2 = 0.0; s.d3 = 0.0; s.d4 = 0.0;
  s.t3cof = 0.0; s.t4cof = 0.0; s.t5cof = 0.0;
  if (s.isimp !== 1) {
    const cc1sq = s.cc1 * s.cc1;
    s.d2 = 4.0 * ao * tsi * cc1sq;
    const tempd = s.d2 * tsi * s.cc1 / 3.0;
    s.d3 = (17.0 * ao + sfour) * tempd;
    s.d4 = 0.5 * tempd * ao * tsi * (221.0 * ao + 31.0 * sfour) * s.cc1;
    s.t3cof = s.d2 + 2.0 * cc1sq;
    s.t4cof = 0.25 * (3.0 * s.d3 + s.cc1 * (12.0 * s.d2 + 10.0 * cc1sq));
    s.t5cof = 0.2 * (3.0 * s.d4 + 12.0 * s.cc1 * s.d3 + 6.0 * s.d2 * s.d2 +
      15.0 * cc1sq * (2.0 * s.d2 + cc1sq));
  }

  // deep space?
  if ((twopi / no) >= 225.0) {
    s.method = "d";
    s.isimp = 1;
    const tc = 0.0;
    const inclm = s.inclo;
    const d = dscom(s.epochjd - 2433281.5, s.ecco, s.argpo, tc, s.inclo, s.nodeo, no);
    // dpper at init
    const op0 = { ep: s.ecco, inclp: s.inclo, nodep: s.nodeo, argpp: s.argpo, mp: s.mo };
    dpper(d.e3, d.ee2, d.peo, d.pgho, d.pho, d.pinco, d.plo, d.se2, d.se3, d.sgh2,
      d.sgh3, d.sgh4, d.sh2, d.sh3, d.si2, d.si3, d.sl2, d.sl3, d.sl4, s.t, d.xgh2,
      d.xgh3, d.xgh4, d.xh2, d.xh3, d.xi2, d.xi3, d.xl2, d.xl3, d.xl4, d.zmol, d.zmos, "y", op0);
    s.ecco = op0.ep; s.inclo = op0.inclp; s.nodeo = op0.nodep; s.argpo = op0.argpp; s.mo = op0.mp;

    // dsinit: pass scalars via args; it writes secular rates + resonance onto d
    const args = {
      no_unkozai: no, ecco: s.ecco, eccsq: eccsq, inclo: inclm,
      gsto: s.gsto, t: 0.0, mo: s.mo, nodeo: s.nodeo, argpo: s.argpo,
      mdot: s.mdot, nodedot: s.nodedot,
    };
    dsinit(d, tc, s.xpidot, args);
    // store deep-space data on sat; dspace works from a fresh state each call
    s.deep = d;
    s.irez = d.irez;
  }

  // build outputs at t=0 sanity
  sgp4(s, 0.0);
  s.init = true;
  return s;
}

/* ---- propagation ------------------------------------------------------- */
function sgp4(s, tsince) {
  const temp4 = 1.5e-12;
  s.t = tsince;
  s.error = 0;

  let mm, argpm, nodem, em, inclm, nm, xincp, ep, argpp, nodep, mp;
  const xmdf = s.mo + s.mdot * s.t;
  const argpdf = s.argpo + s.argpdot * s.t;
  const nodedf = s.nodeo + s.nodedot * s.t;
  argpm = argpdf; mm = xmdf;
  const t2 = s.t * s.t;
  nodem = nodedf + s.nodecf * t2;
  let tempa = 1.0 - s.cc1 * s.t;
  let tempe = s.bstar * s.cc4 * s.t;
  let templ = s.t2cof * t2;

  if (s.isimp !== 1) {
    const delomg = s.omgcof * s.t;
    const delmtemp = 1.0 + s.eta * Math.cos(xmdf);
    const delm = s.xmcof * (delmtemp * delmtemp * delmtemp - s.delmo);
    const temp = delomg + delm;
    mm = xmdf + temp;
    argpm = argpdf - temp;
    const t3 = t2 * s.t, t4 = t3 * s.t;
    tempe = tempe + s.bstar * s.cc5 * (Math.sin(mm) - s.sinmao);
    templ = templ + s.t3cof * t3 + t4 * (s.t4cof + s.t * s.t5cof);
  }

  nm = s.no_unkozai; em = s.ecco; inclm = s.inclo;

  if (s.method === "d") {
    const op = {
      t: s.t, em: s.ecco, inclm: s.inclo, nm: s.no_unkozai,
      argpm: argpm, nodem: nodem, mm: mm,
      no_unkozai: s.no_unkozai, gsto: s.gsto, argpo: s.argpo, argpdot: s.argpdot,
      mo: s.mo, nodeo: s.nodeo,
    };
    dspace(s.deep, s.t, op);
    nm = op.nm; em = op.em; inclm = op.inclm;
    argpm = op.argpm; nodem = op.nodem; mm = op.mm;
  }

  if (nm <= 0.0) { s.error = 2; return { r: null, v: null, error: 2 }; }
  const am = Math.pow(xke / nm, x2o3) * tempa * tempa;
  nm = xke / Math.pow(am, 1.5);
  em = em - tempe;
  if (em >= 1.0 || em < -0.001) { s.error = 1; return { r: null, v: null, error: 1 }; }
  if (em < 1.0e-6) em = 1.0e-6;
  mm = mm + s.no_unkozai * templ;
  let xlm = mm + argpm + nodem;
  const emsq = em * em;
  let temp = 1.0 - emsq;
  nodem = nodem % twopi;
  argpm = argpm % twopi;
  xlm = xlm % twopi;
  mm = (xlm - argpm - nodem) % twopi;

  const sinim = Math.sin(inclm), cosim = Math.cos(inclm);
  ep = em; xincp = inclm; argpp = argpm; nodep = nodem; mp = mm;
  let sinip = sinim, cosip = cosim;

  if (s.method === "d") {
    const op = { ep: ep, inclp: xincp, nodep: nodep, argpp: argpp, mp: mp };
    const d = s.deep;
    dpper(d.e3, d.ee2, d.peo, d.pgho, d.pho, d.pinco, d.plo, d.se2, d.se3, d.sgh2,
      d.sgh3, d.sgh4, d.sh2, d.sh3, d.si2, d.si3, d.sl2, d.sl3, d.sl4, s.t, d.xgh2,
      d.xgh3, d.xgh4, d.xh2, d.xh3, d.xi2, d.xi3, d.xl2, d.xl3, d.xl4, d.zmol, d.zmos, "n", op);
    ep = op.ep; xincp = op.inclp; nodep = op.nodep; argpp = op.argpp; mp = op.mp;
    if (xincp < 0.0) { xincp = -xincp; nodep = nodep + pi; argpp = argpp - pi; }
    if (ep < 0.0 || ep > 1.0) { s.error = 3; return { r: null, v: null, error: 3 }; }
    sinip = Math.sin(xincp); cosip = Math.cos(xincp);
    s.aycof = -0.5 * j3oj2 * sinip;
    if (Math.abs(cosip + 1.0) > 1.5e-12)
      s.xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
    else
      s.xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4;
  }

  const axnl = ep * Math.cos(argpp);
  temp = 1.0 / (am * (1.0 - ep * ep));
  const aynl = ep * Math.sin(argpp) + temp * s.aycof;
  const xl = mp + argpp + nodep + temp * s.xlcof * axnl;

  // Kepler
  let u = (xl - nodep) % twopi;
  let eo1 = u, tem5 = 9999.9, ktr = 1;
  let sineo1 = 0, coseo1 = 0;
  while (Math.abs(tem5) >= 1.0e-12 && ktr <= 10) {
    sineo1 = Math.sin(eo1); coseo1 = Math.cos(eo1);
    tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
    tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
    if (Math.abs(tem5) >= 0.95) tem5 = tem5 > 0.0 ? 0.95 : -0.95;
    eo1 = eo1 + tem5; ktr += 1;
  }

  const ecose = axnl * coseo1 + aynl * sineo1;
  const esine = axnl * sineo1 - aynl * coseo1;
  const el2 = axnl * axnl + aynl * aynl;
  const pl = am * (1.0 - el2);
  if (pl < 0.0) { s.error = 4; return { r: null, v: null, error: 4 }; }

  const rl = am * (1.0 - ecose);
  const rdotl = Math.sqrt(am) * esine / rl;
  const rvdotl = Math.sqrt(pl) / rl;
  const betal = Math.sqrt(1.0 - el2);
  temp = esine / (1.0 + betal);
  const sinu = am / rl * (sineo1 - aynl - axnl * temp);
  const cosu = am / rl * (coseo1 - axnl + aynl * temp);
  let su = Math.atan2(sinu, cosu);
  const sin2u = (cosu + cosu) * sinu;
  const cos2u = 1.0 - 2.0 * sinu * sinu;
  temp = 1.0 / pl;
  const temp1b = 0.5 * j2 * temp;
  const temp2b = temp1b * temp;

  let con41 = s.con41, x1mth2 = s.x1mth2, x7thm1 = s.x7thm1;
  if (s.method === "d") {
    const cosip2 = cosip * cosip;
    con41 = 3.0 * cosip2 - 1.0;
    x1mth2 = 1.0 - cosip2;
    x7thm1 = 7.0 * cosip2 - 1.0;
  }

  let mrt = rl * (1.0 - 1.5 * temp2b * betal * con41) + 0.5 * temp1b * x1mth2 * cos2u;
  su = su - 0.25 * temp2b * x7thm1 * sin2u;
  const xnode = nodep + 1.5 * temp2b * cosip * sin2u;
  const xinc = xincp + 1.5 * temp2b * cosip * sinip * cos2u;
  const mvt = rdotl - nm * temp1b * x1mth2 * sin2u / xke;
  const rvdot = rvdotl + nm * temp1b * (x1mth2 * cos2u + 1.5 * con41) / xke;

  const sinsu = Math.sin(su), cossu = Math.cos(su);
  const snod = Math.sin(xnode), cnod = Math.cos(xnode);
  const sini = Math.sin(xinc), cosi = Math.cos(xinc);
  const xmx = -snod * cosi, xmy = cnod * cosi;
  const ux = xmx * sinsu + cnod * cossu;
  const uy = xmy * sinsu + snod * cossu;
  const uz = sini * sinsu;
  const vx = xmx * cossu - cnod * sinsu;
  const vy = xmy * cossu - snod * sinsu;
  const vz = sini * cossu;

  const r = [mrt * ux * radiusearthkm, mrt * uy * radiusearthkm, mrt * uz * radiusearthkm];
  const vkmpersec = radiusearthkm * xke / 60.0;
  const v = [
    (mvt * ux + rvdot * vx) * vkmpersec,
    (mvt * uy + rvdot * vy) * vkmpersec,
    (mvt * uz + rvdot * vz) * vkmpersec,
  ];
  if (mrt < 1.0) { s.error = 6; return { r, v, error: 6 }; }
  return { r, v, error: 0 };
}

/* ---- helpers ----------------------------------------------------------- */
function jdayFromDate(d) {
  // d: JS Date (UTC). Returns Julian Date.
  return d.getTime() / 86400000.0 + 2440587.5;
}

function sgp4FromOMM(omm) {
  const epochStr = omm.EPOCH.endsWith("Z") ? omm.EPOCH : omm.EPOCH + "Z";
  const epochjd = jdayFromDate(new Date(epochStr));
  return sgp4Init({
    inclo: omm.INCLINATION * deg2rad,
    nodeo: omm.RA_OF_ASC_NODE * deg2rad,
    ecco: omm.ECCENTRICITY,
    argpo: omm.ARG_OF_PERICENTER * deg2rad,
    mo: omm.MEAN_ANOMALY * deg2rad,
    no_kozai: omm.MEAN_MOTION * twopi / 1440.0, // rev/day -> rad/min
    bstar: omm.BSTAR != null ? omm.BSTAR : 0.0,
    epochjd: epochjd,
  });
}

function propagateAt(sat, jsDate) {
  const tsince = (jdayFromDate(jsDate) - sat.epochjd) * 1440.0; // minutes
  return sgp4(sat, tsince);
}

// TEME position -> geodetic lat/lon (deg) + alt (km), using GMST.
function eciToGeodetic(r, jsDate) {
  const jd = jdayFromDate(jsDate);
  const gmst = gstime(jd);
  const x = r[0], y = r[1], z = r[2];
  const a = 6378.137, f = 1.0 / 298.257223563;
  const e2 = 2 * f - f * f;
  let lon = Math.atan2(y, x) - gmst;
  lon = ((lon + pi) % twopi + twopi) % twopi - pi;
  const rxy = Math.sqrt(x * x + y * y);
  let lat = Math.atan2(z, rxy), latOld, C = 0, iter = 0;
  do {
    latOld = lat;
    C = 1.0 / Math.sqrt(1.0 - e2 * Math.sin(lat) * Math.sin(lat));
    lat = Math.atan2(z + a * C * e2 * Math.sin(lat), rxy);
  } while (Math.abs(lat - latOld) > 1e-10 && ++iter < 10);
  const alt = rxy / Math.cos(lat) - a * C;
  return { latDeg: lat / deg2rad, lonDeg: lon / deg2rad, altKm: alt, gmst };
}

window.OLSIM_SGP4 = {
    sgp4Init, sgp4, sgp4FromOMM, propagateAt, eciToGeodetic, gstime,
    jdayFromDate, deg2rad, twopi, radiusearthkm,
  };
/* =========================================================================
 * crossings.js — equator (node) crossing solver on top of sgp4.js
 * -------------------------------------------------------------------------
 * Finds ascending / descending equator crossings of the sub-satellite track
 * by scanning the SGP4 geodetic latitude for sign changes and bisecting to
 * the zero. Returns UTC time + crossing longitude, the two numbers an
 * OSCARLOCATOR reference orbit is plotted from.
 * ========================================================================= */



// Sub-satellite latitude (deg) at a JS Date.
function subLat(sat, date) {
  const { r, error } = propagateAt(sat, date);
  if (error || !r) return null;
  return eciToGeodetic(r, date).latDeg;
}

// Full subpoint at a JS Date.
function subPoint(sat, date) {
  const { r, error } = propagateAt(sat, date);
  if (error || !r) return null;
  return eciToGeodetic(r, date);
}

/*
 * findCrossings(sat, startDate, endDate, opts)
 *   opts.node : "ascending" | "descending" | "both"  (default "ascending")
 *   opts.stepSec : coarse scan step in seconds (default 60)
 * Returns array of { date, lonDeg, latDeg(~0), node } sorted by time.
 */
function findCrossings(sat, startDate, endDate, opts = {}) {
  const node = opts.node || "ascending";
  const stepSec = opts.stepSec || 60;
  const stepMs = stepSec * 1000;
  const t0 = startDate.getTime();
  const t1 = endDate.getTime();
  const out = [];

  let prevT = t0;
  let prevLat = subLat(sat, new Date(prevT));

  for (let t = t0 + stepMs; t <= t1; t += stepMs) {
    const lat = subLat(sat, new Date(t));
    if (prevLat == null || lat == null) { prevT = t; prevLat = lat; continue; }

    const ascending = prevLat < 0 && lat >= 0;   // south -> north
    const descending = prevLat > 0 && lat <= 0;  // north -> south

    let want = false, kind = null;
    if ((node === "ascending" || node === "both") && ascending) { want = true; kind = "ascending"; }
    if ((node === "descending" || node === "both") && descending) { want = true; kind = "descending"; }

    if (want) {
      // bisect for lat == 0 between prevT and t
      let lo = prevT, hi = t, loLat = prevLat, hiLat = lat;
      for (let k = 0; k < 60; k++) {
        const mid = (lo + hi) / 2;
        const ml = subLat(sat, new Date(mid));
        if (ml == null) break;
        // keep the sub-interval that brackets the zero
        if ((loLat < 0 && ml < 0) || (loLat > 0 && ml > 0)) { lo = mid; loLat = ml; }
        else { hi = mid; hiLat = ml; }
        if (hi - lo < 1) break; // 1 ms precision
      }
      const crossMs = (lo + hi) / 2;
      const sp = subPoint(sat, new Date(crossMs));
      if (sp) out.push({ date: new Date(crossMs), lonDeg: sp.lonDeg, latDeg: sp.latDeg, node: kind });
    }
    prevT = t; prevLat = lat;
  }
  return out;
}

/*
 * firstCrossingPerUtcDay(sat, startDate, days, opts)
 *   Builds the classic OSCARLOCATOR reference-orbit table: the FIRST
 *   ascending (and optionally descending) crossing of each UTC day.
 * Returns array of { dateUTC(Date midnight), asc:{date,lonDeg}|null,
 *                    desc:{date,lonDeg}|null }
 */
function firstCrossingPerUtcDay(sat, startDate, days, opts = {}) {
  const wantDesc = !!opts.descending;
  const DAY = 86400000;
  // midnight UTC of startDate
  const d0 = Date.UTC(startDate.getUTCFullYear(), startDate.getUTCMonth(), startDate.getUTCDate());
  const node = wantDesc ? "both" : "ascending";
  // scan whole window once, then bucket by day
  const all = findCrossings(sat, new Date(d0 - DAY), new Date(d0 + (days + 1) * DAY), { node, stepSec: opts.stepSec || 60 });
  const asc = all.filter((c) => c.node === "ascending");
  const desc = all.filter((c) => c.node === "descending");

  const rows = [];
  for (let i = 0; i < days; i++) {
    const dayStart = d0 + i * DAY, dayEnd = dayStart + DAY;
    const a = asc.find((c) => c.date.getTime() >= dayStart && c.date.getTime() < dayEnd);
    const dRow = wantDesc ? desc.find((c) => c.date.getTime() >= dayStart && c.date.getTime() < dayEnd) : null;
    rows.push({
      dateUTC: new Date(dayStart),
      asc: a ? { date: a.date, lonDeg: a.lonDeg } : null,
      desc: dRow ? { date: dRow.date, lonDeg: dRow.lonDeg } : null,
    });
  }
  return rows;
}

window.OLSIM_CROSS = { findCrossings, firstCrossingPerUtcDay, subPoint, subLat };
})();

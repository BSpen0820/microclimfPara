// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include "microclimfCppfunctions.h"
using namespace RcppParallel;

// ============================================================
// Micro1Worker: hourly, static veg, data.frame climate
// ============================================================
struct Micro1Worker : public Worker {
    const RMatrix<double> hgt_m, pai_m, x_m, gsmax_m, lref_m, ltra_m, clump_m, leafd_m, paia_m, leafden_m;
    const RMatrix<double> Smin_m, Smax_m, gref_m, soilb_m, Psie_m, Vq_m, Vm_m, Mc_m, rho_m;
    const RMatrix<double> slope_m, aspect_m, svfa_m, tadd_m;
    const RVector<double> wsa_v, hor_v;
    const RVector<double> tc_v, es_v, ea_v, tdew_v, pk_v, Rsw_v, Rdif_v, lwdown_v, u2_v;
    const RVector<double> soilmp_v, Gp_v, umu_v, kp_v, muGp_v, dtrp_v;
    std::vector<double> zend_r, zenr_r, azid_r;
    std::vector<int>    sindex_r, windex_r;
    std::vector<double> Tgp2_r, Tbp2_r;
    std::vector<bool>   out_r;
    double reqhgt, zref, lat, mxtc, mat;
    int rows, cols, tsteps, ndays, hiy;
    bool complete;
    RVector<double> Tz_w, tleaf_w, relhum_w, soilm_w, uz_w;
    RVector<double> Rdirdown_w, Rdifdown_w, Rlwdown_w, Rswup_w, Rlwup_w;

    Micro1Worker(
        NumericMatrix hgt, NumericMatrix pai, NumericMatrix x, NumericMatrix gsmax,
        NumericMatrix lref, NumericMatrix ltra, NumericMatrix clump, NumericMatrix leafd,
        NumericMatrix paia, NumericMatrix leafden,
        NumericMatrix Smin, NumericMatrix Smax, NumericMatrix gref, NumericMatrix soilb,
        NumericMatrix Psie, NumericMatrix Vq, NumericMatrix Vm, NumericMatrix Mc, NumericMatrix rho,
        NumericMatrix slope, NumericMatrix aspect, NumericMatrix svfa, NumericMatrix tadd,
        NumericVector wsa, NumericVector hor,
        NumericVector tc, NumericVector es, NumericVector ea, NumericVector tdew,
        NumericVector pk, NumericVector Rsw, NumericVector Rdif, NumericVector lwdown, NumericVector u2,
        NumericVector soilmp, NumericVector Gp, NumericVector umu,
        NumericVector kp, NumericVector muGp, NumericVector dtrp,
        std::vector<double> zend, std::vector<double> zenr, std::vector<double> azid,
        std::vector<int> sindex, std::vector<int> windex,
        std::vector<double> Tgp2, std::vector<double> Tbp2, std::vector<bool> out,
        double reqhgt_, double zref_, double lat_, double mxtc_, double mat_,
        int rows_, int cols_, int tsteps_, int ndays_, int hiy_, bool complete_,
        NumericVector Tz, NumericVector tleaf, NumericVector relhum, NumericVector soilm,
        NumericVector uz, NumericVector Rdirdown, NumericVector Rdifdown,
        NumericVector Rlwdown, NumericVector Rswup, NumericVector Rlwup
    ) :
        hgt_m(hgt), pai_m(pai), x_m(x), gsmax_m(gsmax), lref_m(lref), ltra_m(ltra),
        clump_m(clump), leafd_m(leafd), paia_m(paia), leafden_m(leafden),
        Smin_m(Smin), Smax_m(Smax), gref_m(gref), soilb_m(soilb), Psie_m(Psie),
        Vq_m(Vq), Vm_m(Vm), Mc_m(Mc), rho_m(rho),
        slope_m(slope), aspect_m(aspect), svfa_m(svfa), tadd_m(tadd),
        wsa_v(wsa), hor_v(hor),
        tc_v(tc), es_v(es), ea_v(ea), tdew_v(tdew), pk_v(pk),
        Rsw_v(Rsw), Rdif_v(Rdif), lwdown_v(lwdown), u2_v(u2),
        soilmp_v(soilmp), Gp_v(Gp), umu_v(umu), kp_v(kp), muGp_v(muGp), dtrp_v(dtrp),
        zend_r(std::move(zend)), zenr_r(std::move(zenr)), azid_r(std::move(azid)),
        sindex_r(std::move(sindex)), windex_r(std::move(windex)),
        Tgp2_r(std::move(Tgp2)), Tbp2_r(std::move(Tbp2)), out_r(std::move(out)),
        reqhgt(reqhgt_), zref(zref_), lat(lat_), mxtc(mxtc_), mat(mat_),
        rows(rows_), cols(cols_), tsteps(tsteps_), ndays(ndays_), hiy(hiy_), complete(complete_),
        Tz_w(Tz), tleaf_w(tleaf), relhum_w(relhum), soilm_w(soilm), uz_w(uz),
        Rdirdown_w(Rdirdown), Rdifdown_w(Rdifdown), Rlwdown_w(Rlwdown),
        Rswup_w(Rswup), Rlwup_w(Rlwup)
    {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t cell = begin; cell < end; ++cell) {
            int i = static_cast<int>(cell) % rows;
            int j = static_cast<int>(cell) / rows;
            if (std::isnan(hgt_m(i, j))) continue;
            tirstruct tir = twostreamdif(pai_m(i,j), paia_m(i,j), x_m(i,j), lref_m(i,j), ltra_m(i,j), clump_m(i,j), gref_m(i,j));
            stompstruct stomp = stomparamsCpp(hgt_m(i,j), lat, x_m(i,j));
            tiwstruct tiw = windtiCpp(hgt_m(i,j), pai_m(i,j));
            soilpstruct spa;
            spa.Smax = Smax_m(i,j); spa.Smin = Smin_m(i,j); spa.soilb = soilb_m(i,j); spa.psi_e = Psie_m(i,j);
            spa.Vq = Vq_m(i,j); spa.Vm = Vm_m(i,j); spa.Mc = Mc_m(i,j); spa.rho = rho_m(i,j);
            soilstruct sp = soilpfun(Vm_m(i,j), Vq_m(i,j), Mc_m(i,j), rho_m(i,j));
            std::vector<double> Tg(tsteps), DD(tsteps);
            solmodel solp;
            for (int dy = 0; dy < ndays; ++dy) {
                double Rmx = -999.9, tmx = -999.0, tmn = 999.0;
                std::vector<double> surfwet(24), radabs(24), soilmday(24);
                std::vector<double> radCsw(24), radClw(24), Rddown(24), Rbdown(24), radLsw(24), radLpar(24);
                std::vector<double> uf(24), uzday(24), gHa(24);
                for (int hr = 0; hr < 24; ++hr) {
                    int k = dy * 24 + hr;
                    int idx = i + rows * j + cols * rows * k;
                    double si = solarindexCpp(slope_m(i,j), aspect_m(i,j), zend_r[k], azid_r[k], true);
                    if (si < 0.0) si = 0.0;
                    double ws = wsa_v[windex_r[k] * rows * cols + j * rows + i];
                    double ha = hor_v[sindex_r[k] * rows * cols + j * rows + i];
                    double sa = (pi / 2.0) - zenr_r[k];
                    if (ha > std::tan(sa)) si = 0.0;
                    double soild = soildCpp(soilmp_v[k], Smin_m(i,j), Smax_m(i,j), tadd_m(i,j));
                    soilmday[hr] = soild;
                    if (out_r[3]) soilm_w[idx] = soild;
                    solp.zend = zend_r[k]; solp.zenr = zenr_r[k];
                    kstruct kpp = cankCpp(zenr_r[k], x_m(i,j), si);
                    tsdirstruct tspdir = twostreamdirCpp(tir.pait, tir.om, tir.a, tir.gma, tir.J, tir.del, tir.h,
                        gref_m(i,j), kpp.kd, tir.u1, tir.S1, tir.D1, tir.D2);
                    radmodel2 radm = twostreamCpp(pai_m(i,j), clump_m(i,j), gref_m(i,j), svfa_m(i,j), si,
                        tc_v[k], Rsw_v[k], Rdif_v[k], lwdown_v[k], solp, kpp, tspdir, tir);
                    radCsw[hr] = radm.radCsw; radClw[hr] = radm.radClw;
                    Rddown[hr] = radm.Rddown; Rbdown[hr] = radm.Rbdown;
                    radLsw[hr] = radm.radLsw; radLpar[hr] = radm.radLpar;
                    if (out_r[5]) Rdirdown_w[idx] = radm.Rbdown;
                    if (out_r[6]) Rdifdown_w[idx] = radm.Rddown;
                    if (out_r[8]) Rswup_w[idx] = radm.Rdup;
                    double reqhgt2 = (reqhgt < 0.00001) ? 0.00001 : reqhgt;
                    windmodel windm = windCpp(reqhgt2, zref, hgt_m(i,j), pai_m(i,j), u2_v[k], umu_v[k], ws, tiw);
                    uf[hr] = windm.uf; uzday[hr] = windm.uz; gHa[hr] = windm.gHa;
                    if (out_r[4]) uz_w[idx] = windm.uz;
                    soilmodelG0 sG0 = soiltempG0(tc_v[k], es_v[k], ea_v[k], pk_v[k], radm.radGsw, radm.radGlw,
                        tdew_v[k], windm.gHa, soild, mxtc, spa);
                    double Rval = std::abs(sG0.Rnet);
                    if (Rmx < Rval) Rmx = Rval;
                    if (tmx < sG0.Tg) tmx = sG0.Tg;
                    if (tmn > sG0.Tg) tmn = sG0.Tg;
                    surfwet[hr] = sG0.surfwet; radabs[hr] = sG0.radabs;
                }
                double dtr = tmx - tmn;
                for (int hr = 0; hr < 24; ++hr) {
                    int k = dy * 24 + hr;
                    int idx = i + rows * j + cols * rows * k;
                    soilmodel Gvars = soiltemp_hrCpp(tc_v[k], es_v[k], ea_v[k], pk_v[k], radabs[hr], surfwet[hr],
                        tdew_v[k], gHa[hr], soilmday[hr], mxtc, Gp_v[k], dtr, dtrp_v[k], muGp_v[k], kp_v[k], Rmx, sp, spa);
                    Tg[k] = Gvars.Tg; DD[k] = Gvars.DD;
                    if (reqhgt >= 0.0) {
                        radmodel2 rvars;
                        rvars.zend = zend_r[k]; rvars.radCsw = radCsw[hr]; rvars.radClw = radClw[hr];
                        rvars.Rddown = Rddown[hr]; rvars.Rbdown = Rbdown[hr];
                        rvars.radLsw = radLsw[hr]; rvars.radLpar = radLpar[hr];
                        windmodel wvars; wvars.uf = uf[hr]; wvars.uz = uzday[hr]; wvars.gHa = gHa[hr];
                        double reqhgt2 = (reqhgt < 0.00001) ? 0.00001 : reqhgt;
                        abovemodel tv = TVaboveground(reqhgt2, zref, tc_v[k], pk_v[k], ea_v[k], es_v[k], tdew_v[k],
                            Rsw_v[k], Rdif_v[k], lwdown_v[k], soilmday[hr], hgt_m(i,j), pai_m(i,j), paia_m(i,j),
                            x_m(i,j), leafd_m(i,j), leafden_m(i,j), Smin_m(i,j), Smax_m(i,j), Psie_m(i,j),
                            soilb_m(i,j), gsmax_m(i,j), mxtc, stomp, tir, rvars, tiw, wvars, Gvars);
                        if (reqhgt > 0.0) { if (out_r[0]) Tz_w[idx] = tv.Tz; }
                        else              { if (out_r[0]) Tz_w[idx] = Tg[k]; }
                        if (out_r[7]) Rlwdown_w[idx] = tv.lwdn;
                        if (out_r[9]) Rlwup_w[idx] = tv.lwup;
                        if (reqhgt > 0.0) {
                            if (out_r[1]) tleaf_w[idx] = tv.tleaf;
                            if (out_r[2]) relhum_w[idx] = tv.rh;
                        }
                    }
                }
            }
            if (reqhgt < 0.0 && out_r[0]) {
                double sumD = 0.0;
                for (int k = 0; k < tsteps; ++k) sumD += DD[k];
                double meanD = sumD / static_cast<double>(tsteps);
                std::vector<double> Tzv = Tbelowgroundv(reqhgt, Tg, Tgp2_r, Tbp2_r, meanD, mat, hiy, complete);
                for (int k = 0; k < tsteps; ++k)
                    Tz_w[i + rows * j + cols * rows * k] = Tzv[k];
            }
        }
    }
};

// [[Rcpp::export]]
List runmicro1Par(DataFrame obstime, DataFrame climdata, DataFrame pointm, List vegp, List soilc,
    double reqhgt, double zref, double lat, double lon, double Sminp, double Smaxp, double tfact,
    bool complete, double mat, std::vector<bool> out, int ncores = 2)
{
    IntegerVector year = obstime["year"]; IntegerVector month = obstime["month"];
    IntegerVector day = obstime["day"];   NumericVector hour = obstime["hour"];
    NumericVector tc = climdata["temp"]; NumericVector es = climdata["es"];
    NumericVector ea = climdata["ea"];   NumericVector tdew = climdata["tdew"];
    NumericVector pk = climdata["pres"]; NumericVector Rsw = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"]; NumericVector lwdown = climdata["lwdown"];
    NumericVector u2 = climdata["windspeed"]; NumericVector wdir = climdata["winddir"];
    NumericVector soilmp = pointm["soilm"]; NumericVector Tgp = pointm["Tg"];
    NumericVector Tbp = pointm["Tbp"];      NumericVector Gp = pointm["G"];
    NumericVector umu = pointm["umu"]; NumericVector kp = pointm["kp"];
    NumericVector muGp = pointm["muGp"]; NumericVector dtrp = pointm["dtrp"];
    NumericMatrix hgt = vegp["hgt"]; NumericMatrix pai = vegp["pai"];
    NumericMatrix x = vegp["x"]; NumericMatrix gsmax = vegp["gsmax"];
    NumericMatrix lref = vegp["leafr"]; NumericMatrix ltra = vegp["leaft"];
    NumericMatrix clump = vegp["clump"]; NumericMatrix leafd = vegp["leafd"];
    NumericMatrix paia = vegp["paia"]; NumericMatrix leafden = vegp["leafden"];
    NumericMatrix Smin = soilc["Smin"]; NumericMatrix Smax = soilc["Smax"];
    NumericMatrix gref = soilc["gref"]; NumericMatrix soilb = soilc["soilb"];
    NumericMatrix Psie = soilc["Psie"]; NumericMatrix Vq = soilc["Vq"];
    NumericMatrix Vm = soilc["Vm"]; NumericMatrix Mc = soilc["Mc"]; NumericMatrix rho = soilc["rho"];
    NumericMatrix slope = soilc["slope"]; NumericMatrix aspect = soilc["aspect"];
    NumericMatrix twi = soilc["twi"]; NumericMatrix svfa = soilc["svfa"];
    NumericVector wsa = soilc["wsa"]; NumericVector hor = soilc["hor"];
    int rows = hgt.nrow(); int cols = hgt.ncol();
    int tsteps = tc.size(); int ndays = tsteps / 24;
    IntegerVector dim = {rows, cols, tsteps}; int n = rows * cols * tsteps;
    std::vector<double> zend_v(tsteps), zenr_v(tsteps), azid_v(tsteps);
    std::vector<int> sindex_v(tsteps), windex_v(tsteps);
    double mxtc = -273.15;
    for (int k = 0; k < tsteps; ++k) {
        solmodel sp = solpositionCpp(lat, lon, year[k], month[k], day[k], hour[k]);
        zend_v[k] = sp.zend; zenr_v[k] = sp.zenr; azid_v[k] = sp.azid;
        sindex_v[k] = static_cast<int>(std::round(sp.azid / 15.0)) % 24;
        windex_v[k] = static_cast<int>(std::round(wdir[k] / 45)) % 8;
        if (tc[k] > mxtc) mxtc = tc[k];
    }
    int hiy = (year[0] % 4 == 0) ? 366 * 24 : 365 * 24;
    NumericMatrix tadd = soildCppm(twi, Sminp, Smaxp, tfact);
    std::vector<double> Tgp2 = as<std::vector<double>>(Tgp);
    std::vector<double> Tbp2 = as<std::vector<double>>(Tbp);
    // Only allocate the output arrays actually requested in `out`. All worker
    // writes are guarded by out_r[k], so unrequested outputs (left empty) are
    // never indexed. Mirrors serial runmicro1-4Cpp; cuts peak output memory by
    // up to 10x when few layers are requested.
    NumericVector Tz, tleaf, relhum, soilmv, uz, Rdirdown, Rdifdown, Rlwdown, Rswup, Rlwup;
    if (out[0]) Tz       = NumericVector(n, NA_REAL);
    if (out[1]) tleaf    = NumericVector(n, NA_REAL);
    if (out[2]) relhum   = NumericVector(n, NA_REAL);
    if (out[3]) soilmv   = NumericVector(n, NA_REAL);
    if (out[4]) uz       = NumericVector(n, NA_REAL);
    if (out[5]) Rdirdown = NumericVector(n, NA_REAL);
    if (out[6]) Rdifdown = NumericVector(n, NA_REAL);
    if (out[7]) Rlwdown  = NumericVector(n, NA_REAL);
    if (out[8]) Rswup    = NumericVector(n, NA_REAL);
    if (out[9]) Rlwup    = NumericVector(n, NA_REAL);
    Micro1Worker worker(
        hgt, pai, x, gsmax, lref, ltra, clump, leafd, paia, leafden,
        Smin, Smax, gref, soilb, Psie, Vq, Vm, Mc, rho, slope, aspect, svfa, tadd,
        wsa, hor, tc, es, ea, tdew, pk, Rsw, Rdif, lwdown, u2,
        soilmp, Gp, umu, kp, muGp, dtrp,
        zend_v, zenr_v, azid_v, sindex_v, windex_v, Tgp2, Tbp2, out,
        reqhgt, zref, lat, mxtc, mat, rows, cols, tsteps, ndays, hiy, complete,
        Tz, tleaf, relhum, soilmv, uz, Rdirdown, Rdifdown, Rlwdown, Rswup, Rlwup);
    parallelFor(0, rows * cols, worker, 1, ncores);
    if (out[0]) Tz.attr("dim") = dim;
    if (out[1]) tleaf.attr("dim") = dim;
    if (out[2]) relhum.attr("dim") = dim;
    if (out[3]) soilmv.attr("dim") = dim;
    if (out[4]) uz.attr("dim") = dim;
    if (out[5]) Rdirdown.attr("dim") = dim;
    if (out[6]) Rdifdown.attr("dim") = dim;
    if (out[7]) Rlwdown.attr("dim") = dim;
    if (out[8]) Rswup.attr("dim") = dim;
    if (out[9]) Rlwup.attr("dim") = dim;
    Rcpp::List outp;
    if (out[0]) outp["Tz"] = Tz;
    if (out[1]) outp["tleaf"] = tleaf;
    if (out[2]) outp["relhum"] = relhum;
    if (out[3]) outp["soilm"] = soilmv;
    if (out[4]) outp["windspeed"] = uz;
    if (out[5]) outp["Rdirdown"] = Rdirdown;
    if (out[6]) outp["Rdifdown"] = Rdifdown;
    if (out[7]) outp["Rlwdown"] = Rlwdown;
    if (out[8]) outp["Rswup"] = Rswup;
    if (out[9]) outp["Rlwup"] = Rlwup;
    return outp;
}

// ============================================================
// Micro2Worker: hourly, static veg, array climate
// ============================================================
struct Micro2Worker : public Worker {
    const RMatrix<double> hgt_m, pai_m, x_m, gsmax_m, lref_m, ltra_m, clump_m, leafd_m, paia_m, leafden_m;
    const RMatrix<double> Smin_m, Smax_m, gref_m, soilb_m, Psie_m, Vq_m, Vm_m, Mc_m, rho_m;
    const RMatrix<double> slope_m, aspect_m, svfa_m, tadd_m, lats_m, lons_m;
    const RVector<double> wsa_v, hor_v;
    const RVector<double> tc_v, es_v, ea_v, tdew_v, pk_v, Rsw_v, Rdif_v, lwdown_v, u2_v;
    const RVector<double> soilmp_v, Tgp_v, Tbp_v, Gp_v, umu_v, kp_v, muGp_v, dtrp_v;
    std::vector<int>  windex_r;
    std::vector<int>  year_r, month_r, day_r;
    std::vector<double> hour_r;
    std::vector<bool> out_r;
    double reqhgt, zref, mat;
    int rows, cols, tsteps, ndays, hiy;
    bool complete;
    RVector<double> Tz_w, tleaf_w, relhum_w, soilm_w, uz_w;
    RVector<double> Rdirdown_w, Rdifdown_w, Rlwdown_w, Rswup_w, Rlwup_w;

    Micro2Worker(
        NumericMatrix hgt, NumericMatrix pai, NumericMatrix x, NumericMatrix gsmax,
        NumericMatrix lref, NumericMatrix ltra, NumericMatrix clump, NumericMatrix leafd,
        NumericMatrix paia, NumericMatrix leafden,
        NumericMatrix Smin, NumericMatrix Smax, NumericMatrix gref, NumericMatrix soilb,
        NumericMatrix Psie, NumericMatrix Vq, NumericMatrix Vm, NumericMatrix Mc, NumericMatrix rho,
        NumericMatrix slope, NumericMatrix aspect, NumericMatrix svfa, NumericMatrix tadd,
        NumericMatrix lats, NumericMatrix lons,
        NumericVector wsa, NumericVector hor,
        NumericVector tc, NumericVector es, NumericVector ea, NumericVector tdew,
        NumericVector pk, NumericVector Rsw, NumericVector Rdif, NumericVector lwdown, NumericVector u2,
        NumericVector soilmp, NumericVector Tgp, NumericVector Tbp, NumericVector Gp,
        NumericVector umu, NumericVector kp, NumericVector muGp, NumericVector dtrp,
        std::vector<int> windex, std::vector<int> year, std::vector<int> month,
        std::vector<int> day, std::vector<double> hour,
        std::vector<bool> out,
        double reqhgt_, double zref_, double mat_,
        int rows_, int cols_, int tsteps_, int ndays_, int hiy_, bool complete_,
        NumericVector Tz, NumericVector tleaf, NumericVector relhum, NumericVector soilm,
        NumericVector uz, NumericVector Rdirdown, NumericVector Rdifdown,
        NumericVector Rlwdown, NumericVector Rswup, NumericVector Rlwup
    ) :
        hgt_m(hgt), pai_m(pai), x_m(x), gsmax_m(gsmax), lref_m(lref), ltra_m(ltra),
        clump_m(clump), leafd_m(leafd), paia_m(paia), leafden_m(leafden),
        Smin_m(Smin), Smax_m(Smax), gref_m(gref), soilb_m(soilb), Psie_m(Psie),
        Vq_m(Vq), Vm_m(Vm), Mc_m(Mc), rho_m(rho),
        slope_m(slope), aspect_m(aspect), svfa_m(svfa), tadd_m(tadd), lats_m(lats), lons_m(lons),
        wsa_v(wsa), hor_v(hor),
        tc_v(tc), es_v(es), ea_v(ea), tdew_v(tdew), pk_v(pk),
        Rsw_v(Rsw), Rdif_v(Rdif), lwdown_v(lwdown), u2_v(u2),
        soilmp_v(soilmp), Tgp_v(Tgp), Tbp_v(Tbp), Gp_v(Gp),
        umu_v(umu), kp_v(kp), muGp_v(muGp), dtrp_v(dtrp),
        windex_r(std::move(windex)), year_r(std::move(year)), month_r(std::move(month)),
        day_r(std::move(day)), hour_r(std::move(hour)), out_r(std::move(out)),
        reqhgt(reqhgt_), zref(zref_), mat(mat_),
        rows(rows_), cols(cols_), tsteps(tsteps_), ndays(ndays_), hiy(hiy_), complete(complete_),
        Tz_w(Tz), tleaf_w(tleaf), relhum_w(relhum), soilm_w(soilm), uz_w(uz),
        Rdirdown_w(Rdirdown), Rdifdown_w(Rdifdown), Rlwdown_w(Rlwdown),
        Rswup_w(Rswup), Rlwup_w(Rlwup)
    {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t cell = begin; cell < end; ++cell) {
            int i = static_cast<int>(cell) % rows;
            int j = static_cast<int>(cell) / rows;
            if (std::isnan(hgt_m(i, j))) continue;
            tirstruct tir = twostreamdif(pai_m(i,j), paia_m(i,j), x_m(i,j), lref_m(i,j), ltra_m(i,j), clump_m(i,j), gref_m(i,j));
            stompstruct stomp = stomparamsCpp(hgt_m(i,j), lats_m(i,j), x_m(i,j));
            tiwstruct tiw = windtiCpp(hgt_m(i,j), pai_m(i,j));
            soilpstruct spa;
            spa.Smax = Smax_m(i,j); spa.Smin = Smin_m(i,j); spa.soilb = soilb_m(i,j); spa.psi_e = Psie_m(i,j);
            spa.Vq = Vq_m(i,j); spa.Vm = Vm_m(i,j); spa.Mc = Mc_m(i,j); spa.rho = rho_m(i,j);
            soilstruct sp = soilpfun(Vm_m(i,j), Vq_m(i,j), Mc_m(i,j), rho_m(i,j));
            double mxtc = -273.15;
            for (int k = 0; k < tsteps; ++k) {
                double v = tc_v[i + rows * j + cols * rows * k];
                if (v > mxtc) mxtc = v;
            }
            std::vector<double> Tg(tsteps), DD(tsteps);
            for (int dy = 0; dy < ndays; ++dy) {
                double Rmx = -999.9, tmx = -999.0, tmn = 999.0;
                std::vector<double> surfwet(24), radabs(24), soilmday(24);
                std::vector<double> radCsw(24), radClw(24), Rddown(24), Rbdown(24), radLsw(24), radLpar(24);
                std::vector<double> uf(24), uzday(24), gHa(24), zendday(24);
                for (int hr = 0; hr < 24; ++hr) {
                    int k = dy * 24 + hr;
                    int idx = i + rows * j + cols * rows * k;
                    solmodel solp = solpositionCpp(lats_m(i,j), lons_m(i,j), year_r[k], month_r[k], day_r[k], hour_r[k]);
                    zendday[hr] = solp.zend;
                    double si = solarindexCpp(slope_m(i,j), aspect_m(i,j), solp.zend, solp.azid);
                    if (si < 0.0) si = 0.0;
                    int sindx = static_cast<int>(std::round(solp.azid / 15)) % 24;
                    double ha = hor_v[sindx * rows * cols + j * rows + i];
                    double sa = (pi / 2.0) - solp.zenr;
                    if (ha > std::tan(sa)) si = 0.0;
                    double ws = wsa_v[windex_r[k] * rows * cols + j * rows + i];
                    double soild = soildCpp(soilmp_v[idx], Smin_m(i,j), Smax_m(i,j), tadd_m(i,j));
                    soilmday[hr] = soild;
                    if (out_r[3]) soilm_w[idx] = soild;
                    kstruct kpp = cankCpp(solp.zenr, x_m(i,j), si);
                    tsdirstruct tspdir = twostreamdirCpp(tir.pait, tir.om, tir.a, tir.gma, tir.J, tir.del, tir.h,
                        gref_m(i,j), kpp.kd, tir.u1, tir.S1, tir.D1, tir.D2);
                    radmodel2 radm = twostreamCpp(pai_m(i,j), clump_m(i,j), gref_m(i,j), svfa_m(i,j), si,
                        tc_v[idx], Rsw_v[idx], Rdif_v[idx], lwdown_v[idx], solp, kpp, tspdir, tir);
                    radCsw[hr] = radm.radCsw; radClw[hr] = radm.radClw;
                    Rddown[hr] = radm.Rddown; Rbdown[hr] = radm.Rbdown;
                    radLsw[hr] = radm.radLsw; radLpar[hr] = radm.radLpar;
                    if (out_r[5]) Rdirdown_w[idx] = radm.Rbdown;
                    if (out_r[6]) Rdifdown_w[idx] = radm.Rddown;
                    if (out_r[8]) Rswup_w[idx] = radm.Rdup;
                    double reqhgt2 = (reqhgt < 0.00001) ? 0.00001 : reqhgt;
                    windmodel windm = windCpp(reqhgt2, zref, hgt_m(i,j), pai_m(i,j), u2_v[idx], umu_v[idx], ws, tiw);
                    uf[hr] = windm.uf; uzday[hr] = windm.uz; gHa[hr] = windm.gHa;
                    if (out_r[4]) uz_w[idx] = windm.uz;
                    soilmodelG0 sG0 = soiltempG0(tc_v[idx], es_v[idx], ea_v[idx], pk_v[idx], radm.radGsw, radm.radGlw,
                        tdew_v[idx], windm.gHa, soild, mxtc, spa);
                    double Rval = std::abs(sG0.Rnet);
                    if (Rmx < Rval) Rmx = Rval;
                    if (tmx < sG0.Tg) tmx = sG0.Tg;
                    if (tmn > sG0.Tg) tmn = sG0.Tg;
                    surfwet[hr] = sG0.surfwet; radabs[hr] = sG0.radabs;
                }
                double dtr = tmx - tmn;
                for (int hr = 0; hr < 24; ++hr) {
                    int k = dy * 24 + hr;
                    int idx = i + rows * j + cols * rows * k;
                    soilmodel Gvars = soiltemp_hrCpp(tc_v[idx], es_v[idx], ea_v[idx], pk_v[idx], radabs[hr], surfwet[hr],
                        tdew_v[idx], gHa[hr], soilmday[hr], mxtc, Gp_v[idx], dtr, dtrp_v[idx], muGp_v[idx], kp_v[idx], Rmx, sp, spa);
                    Tg[k] = Gvars.Tg; DD[k] = Gvars.DD;
                    if (reqhgt >= 0.0) {
                        radmodel2 rvars;
                        rvars.zend = zendday[hr]; rvars.radCsw = radCsw[hr]; rvars.radClw = radClw[hr];
                        rvars.Rddown = Rddown[hr]; rvars.Rbdown = Rbdown[hr];
                        rvars.radLsw = radLsw[hr]; rvars.radLpar = radLpar[hr];
                        windmodel wvars; wvars.uf = uf[hr]; wvars.uz = uzday[hr]; wvars.gHa = gHa[hr];
                        double reqhgt2 = (reqhgt < 0.00001) ? 0.00001 : reqhgt;
                        abovemodel tv = TVaboveground(reqhgt2, zref, tc_v[idx], pk_v[idx], ea_v[idx], es_v[idx], tdew_v[idx],
                            Rsw_v[idx], Rdif_v[idx], lwdown_v[idx], soilmday[hr], hgt_m(i,j), pai_m(i,j), paia_m(i,j),
                            x_m(i,j), leafd_m(i,j), leafden_m(i,j), Smin_m(i,j), Smax_m(i,j), Psie_m(i,j),
                            soilb_m(i,j), gsmax_m(i,j), mxtc, stomp, tir, rvars, tiw, wvars, Gvars);
                        if (reqhgt > 0.0) { if (out_r[0]) Tz_w[idx] = tv.Tz; }
                        else              { if (out_r[0]) Tz_w[idx] = Tg[k]; }
                        if (out_r[7]) Rlwdown_w[idx] = tv.lwdn;
                        if (out_r[9]) Rlwup_w[idx] = tv.lwup;
                        if (reqhgt > 0.0) {
                            if (out_r[1]) tleaf_w[idx] = tv.tleaf;
                            if (out_r[2]) relhum_w[idx] = tv.rh;
                        }
                    }
                }
            }
            if (reqhgt < 0.0 && out_r[0]) {
                std::vector<double> Tgpv(tsteps), Tbpv(tsteps);
                double sumD = 0.0;
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k;
                    Tgpv[k] = Tgp_v[idx]; Tbpv[k] = Tbp_v[idx];
                    sumD += DD[k];
                }
                double meanD = sumD / static_cast<double>(tsteps);
                std::vector<double> Tzv = Tbelowgroundv(reqhgt, Tg, Tgpv, Tbpv, meanD, mat, hiy, complete);
                for (int k = 0; k < tsteps; ++k)
                    Tz_w[i + rows * j + cols * rows * k] = Tzv[k];
            }
        }
    }
};

// [[Rcpp::export]]
List runmicro2Par(DataFrame obstime, List climdata, List pointm, List vegp, List soilc,
    double reqhgt, double zref, NumericMatrix lats, NumericMatrix lons,
    double Sminp, double Smaxp, double tfact, bool complete, double mat,
    std::vector<bool> out, int ncores = 2)
{
    IntegerVector year = obstime["year"]; IntegerVector month = obstime["month"];
    IntegerVector day = obstime["day"];   NumericVector hour = obstime["hour"];
    NumericVector tc = climdata["tc"]; NumericVector es = climdata["es"];
    NumericVector ea = climdata["ea"]; NumericVector tdew = climdata["tdew"];
    NumericVector pk = climdata["pk"]; NumericVector Rsw = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"]; NumericVector lwdown = climdata["lwdown"];
    NumericVector u2 = climdata["windspeed"]; NumericVector wdir = climdata["winddir"];
    NumericVector soilmp = pointm["soilm"]; NumericVector Tgp = pointm["Tg"];
    NumericVector Tbp = pointm["Tbp"]; NumericVector Gp = pointm["Gp"];
    NumericVector umu = pointm["umu"]; NumericVector kp = pointm["kp"];
    NumericVector muGp = pointm["muGp"]; NumericVector dtrp = pointm["dtrp"];
    NumericMatrix hgt = vegp["hgt"]; NumericMatrix pai = vegp["pai"];
    NumericMatrix x = vegp["x"]; NumericMatrix gsmax = vegp["gsmax"];
    NumericMatrix lref = vegp["leafr"]; NumericMatrix ltra = vegp["leaft"];
    NumericMatrix clump = vegp["clump"]; NumericMatrix leafd = vegp["leafd"];
    NumericMatrix paia = vegp["paia"]; NumericMatrix leafden = vegp["leafden"];
    NumericMatrix Smin = soilc["Smin"]; NumericMatrix Smax = soilc["Smax"];
    NumericMatrix gref = soilc["gref"]; NumericMatrix soilb = soilc["soilb"];
    NumericMatrix Psie = soilc["Psie"]; NumericMatrix Vq = soilc["Vq"];
    NumericMatrix Vm = soilc["Vm"]; NumericMatrix Mc = soilc["Mc"]; NumericMatrix rho = soilc["rho"];
    NumericMatrix slope = soilc["slope"]; NumericMatrix aspect = soilc["aspect"];
    NumericMatrix twi = soilc["twi"]; NumericMatrix svfa = soilc["svfa"];
    NumericVector wsa = soilc["wsa"]; NumericVector hor = soilc["hor"];
    int rows = hgt.nrow(); int cols = hgt.ncol();
    int tsteps = year.size(); int ndays = tsteps / 24;
    IntegerVector dim = {rows, cols, tsteps}; int n = rows * cols * tsteps;
    std::vector<int> windex_v(tsteps);
    for (int k = 0; k < tsteps; ++k)
        windex_v[k] = static_cast<int>(std::round(wdir[k] / 45)) % 8;
    int hiy = (year[0] % 4 == 0) ? 366 * 24 : 365 * 24;
    NumericMatrix tadd = soildCppm(twi, Sminp, Smaxp, tfact);
    std::vector<int> yv(tsteps), mv(tsteps), dv(tsteps);
    std::vector<double> hv(tsteps);
    for (int k = 0; k < tsteps; ++k) { yv[k]=year[k]; mv[k]=month[k]; dv[k]=day[k]; hv[k]=hour[k]; }
    // Only allocate the output arrays actually requested in `out`. All worker
    // writes are guarded by out_r[k], so unrequested outputs (left empty) are
    // never indexed. Mirrors serial runmicro1-4Cpp; cuts peak output memory by
    // up to 10x when few layers are requested.
    NumericVector Tz, tleaf, relhum, soilmv, uz, Rdirdown, Rdifdown, Rlwdown, Rswup, Rlwup;
    if (out[0]) Tz       = NumericVector(n, NA_REAL);
    if (out[1]) tleaf    = NumericVector(n, NA_REAL);
    if (out[2]) relhum   = NumericVector(n, NA_REAL);
    if (out[3]) soilmv   = NumericVector(n, NA_REAL);
    if (out[4]) uz       = NumericVector(n, NA_REAL);
    if (out[5]) Rdirdown = NumericVector(n, NA_REAL);
    if (out[6]) Rdifdown = NumericVector(n, NA_REAL);
    if (out[7]) Rlwdown  = NumericVector(n, NA_REAL);
    if (out[8]) Rswup    = NumericVector(n, NA_REAL);
    if (out[9]) Rlwup    = NumericVector(n, NA_REAL);
    Micro2Worker worker(
        hgt, pai, x, gsmax, lref, ltra, clump, leafd, paia, leafden,
        Smin, Smax, gref, soilb, Psie, Vq, Vm, Mc, rho, slope, aspect, svfa, tadd, lats, lons,
        wsa, hor, tc, es, ea, tdew, pk, Rsw, Rdif, lwdown, u2,
        soilmp, Tgp, Tbp, Gp, umu, kp, muGp, dtrp,
        windex_v, yv, mv, dv, hv, out,
        reqhgt, zref, mat, rows, cols, tsteps, ndays, hiy, complete,
        Tz, tleaf, relhum, soilmv, uz, Rdirdown, Rdifdown, Rlwdown, Rswup, Rlwup);
    parallelFor(0, rows * cols, worker, 1, ncores);
    if (out[0]) Tz.attr("dim") = dim;
    if (out[1]) tleaf.attr("dim") = dim;
    if (out[2]) relhum.attr("dim") = dim;
    if (out[3]) soilmv.attr("dim") = dim;
    if (out[4]) uz.attr("dim") = dim;
    if (out[5]) Rdirdown.attr("dim") = dim;
    if (out[6]) Rdifdown.attr("dim") = dim;
    if (out[7]) Rlwdown.attr("dim") = dim;
    if (out[8]) Rswup.attr("dim") = dim;
    if (out[9]) Rlwup.attr("dim") = dim;
    Rcpp::List outp;
    if (out[0]) outp["Tz"] = Tz;
    if (out[1]) outp["tleaf"] = tleaf;
    if (out[2]) outp["relhum"] = relhum;
    if (out[3]) outp["soilm"] = soilmv;
    if (out[4]) outp["windspeed"] = uz;
    if (out[5]) outp["Rdirdown"] = Rdirdown;
    if (out[6]) outp["Rdifdown"] = Rdifdown;
    if (out[7]) outp["Rlwdown"] = Rlwdown;
    if (out[8]) outp["Rswup"] = Rswup;
    if (out[9]) outp["Rlwup"] = Rlwup;
    return outp;
}

// ============================================================
// Micro3Worker: hourly, time-variant veg, data.frame climate
// ============================================================
struct Micro3Worker : public Worker {
    const RVector<double> hgt_v, pai_v, x_v, gsmax_v, lref_v, ltra_v, clump_v, leafd_v, paia_v, leafden_v;
    const RMatrix<double> Smin_m, Smax_m, gref_m, soilb_m, Psie_m, Vq_m, Vm_m, Mc_m, rho_m;
    const RMatrix<double> slope_m, aspect_m, svfa_m, tadd_m;
    const RVector<double> wsa_v, hor_v;
    const RVector<double> tc_v, es_v, ea_v, tdew_v, pk_v, Rsw_v, Rdif_v, lwdown_v, u2_v;
    const RVector<double> soilmp_v, Gp_v, umu_v, kp_v, muGp_v, dtrp_v;
    std::vector<double> zend_r, zenr_r, azid_r;
    std::vector<int>    sindex_r, windex_r;
    std::vector<double> Tgp2_r, Tbp2_r;
    std::vector<bool>   out_r;
    std::vector<int>    lyr_r, st_r, ndays_r;
    int nlyrs;
    double reqhgt, zref, lat, mxtc, mat;
    int rows, cols, tsteps, hiy;
    bool complete;
    RVector<double> Tz_w, tleaf_w, relhum_w, soilm_w, uz_w;
    RVector<double> Rdirdown_w, Rdifdown_w, Rlwdown_w, Rswup_w, Rlwup_w;

    Micro3Worker(
        NumericVector hgt, NumericVector pai, NumericVector x, NumericVector gsmax,
        NumericVector lref, NumericVector ltra, NumericVector clump, NumericVector leafd,
        NumericVector paia, NumericVector leafden,
        NumericMatrix Smin, NumericMatrix Smax, NumericMatrix gref, NumericMatrix soilb,
        NumericMatrix Psie, NumericMatrix Vq, NumericMatrix Vm, NumericMatrix Mc, NumericMatrix rho,
        NumericMatrix slope, NumericMatrix aspect, NumericMatrix svfa, NumericMatrix tadd,
        NumericVector wsa, NumericVector hor,
        NumericVector tc, NumericVector es, NumericVector ea, NumericVector tdew,
        NumericVector pk, NumericVector Rsw, NumericVector Rdif, NumericVector lwdown, NumericVector u2,
        NumericVector soilmp, NumericVector Gp, NumericVector umu,
        NumericVector kp, NumericVector muGp, NumericVector dtrp,
        std::vector<double> zend, std::vector<double> zenr, std::vector<double> azid,
        std::vector<int> sindex, std::vector<int> windex,
        std::vector<double> Tgp2, std::vector<double> Tbp2, std::vector<bool> out,
        std::vector<int> lyr, std::vector<int> st, std::vector<int> ndays, int nlyrs_,
        double reqhgt_, double zref_, double lat_, double mxtc_, double mat_,
        int rows_, int cols_, int tsteps_, int hiy_, bool complete_,
        NumericVector Tz, NumericVector tleaf, NumericVector relhum, NumericVector soilm,
        NumericVector uz, NumericVector Rdirdown, NumericVector Rdifdown,
        NumericVector Rlwdown, NumericVector Rswup, NumericVector Rlwup
    ) :
        hgt_v(hgt), pai_v(pai), x_v(x), gsmax_v(gsmax), lref_v(lref), ltra_v(ltra),
        clump_v(clump), leafd_v(leafd), paia_v(paia), leafden_v(leafden),
        Smin_m(Smin), Smax_m(Smax), gref_m(gref), soilb_m(soilb), Psie_m(Psie),
        Vq_m(Vq), Vm_m(Vm), Mc_m(Mc), rho_m(rho),
        slope_m(slope), aspect_m(aspect), svfa_m(svfa), tadd_m(tadd),
        wsa_v(wsa), hor_v(hor),
        tc_v(tc), es_v(es), ea_v(ea), tdew_v(tdew), pk_v(pk),
        Rsw_v(Rsw), Rdif_v(Rdif), lwdown_v(lwdown), u2_v(u2),
        soilmp_v(soilmp), Gp_v(Gp), umu_v(umu), kp_v(kp), muGp_v(muGp), dtrp_v(dtrp),
        zend_r(std::move(zend)), zenr_r(std::move(zenr)), azid_r(std::move(azid)),
        sindex_r(std::move(sindex)), windex_r(std::move(windex)),
        Tgp2_r(std::move(Tgp2)), Tbp2_r(std::move(Tbp2)), out_r(std::move(out)),
        lyr_r(std::move(lyr)), st_r(std::move(st)), ndays_r(std::move(ndays)), nlyrs(nlyrs_),
        reqhgt(reqhgt_), zref(zref_), lat(lat_), mxtc(mxtc_), mat(mat_),
        rows(rows_), cols(cols_), tsteps(tsteps_), hiy(hiy_), complete(complete_),
        Tz_w(Tz), tleaf_w(tleaf), relhum_w(relhum), soilm_w(soilm), uz_w(uz),
        Rdirdown_w(Rdirdown), Rdifdown_w(Rdifdown), Rlwdown_w(Rlwdown),
        Rswup_w(Rswup), Rlwup_w(Rlwup)
    {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t cell = begin; cell < end; ++cell) {
            int i = static_cast<int>(cell) % rows;
            int j = static_cast<int>(cell) / rows;
            int base = i + rows * j;
            if (std::isnan(hgt_v[base])) continue;
            soilpstruct spa;
            spa.Smax = Smax_m(i,j); spa.Smin = Smin_m(i,j); spa.soilb = soilb_m(i,j); spa.psi_e = Psie_m(i,j);
            spa.Vq = Vq_m(i,j); spa.Vm = Vm_m(i,j); spa.Mc = Mc_m(i,j); spa.rho = rho_m(i,j);
            soilstruct sp = soilpfun(Vm_m(i,j), Vq_m(i,j), Mc_m(i,j), rho_m(i,j));
            std::vector<double> Tg(tsteps), DD(tsteps);
            solmodel solp;
            for (int lyr = 0; lyr < nlyrs; ++lyr) {
                int idxl = base + cols * rows * lyr;
                tirstruct tir = twostreamdif(pai_v[idxl], paia_v[idxl], x_v[idxl], lref_v[idxl], ltra_v[idxl], clump_v[idxl], gref_m(i,j));
                stompstruct stomp = stomparamsCpp(hgt_v[idxl], lat, x_v[idxl]);
                tiwstruct tiw = windtiCpp(hgt_v[idxl], pai_v[idxl]);
                for (int dy = 0; dy < ndays_r[lyr]; ++dy) {
                    double Rmx = -999.9, tmx = -999.0, tmn = 999.0;
                    std::vector<double> surfwet(24), radabs(24), soilmday(24);
                    std::vector<double> radCsw(24), radClw(24), Rddown(24), Rbdown(24), radLsw(24), radLpar(24);
                    std::vector<double> uf(24), uzday(24), gHa(24);
                    for (int hr = 0; hr < 24; ++hr) {
                        int k = dy * 24 + hr + st_r[lyr];
                        int idx = base + cols * rows * k;
                        double si = solarindexCpp(slope_m(i,j), aspect_m(i,j), zend_r[k], azid_r[k], true);
                        if (si < 0.0) si = 0.0;
                        double ws = wsa_v[windex_r[k] * rows * cols + j * rows + i];
                        double ha = hor_v[sindex_r[k] * rows * cols + j * rows + i];
                        double sa = (pi / 2.0) - zenr_r[k];
                        if (ha > std::tan(sa)) si = 0.0;
                        double soild = soildCpp(soilmp_v[k], Smin_m(i,j), Smax_m(i,j), tadd_m(i,j));
                        soilmday[hr] = soild;
                        if (out_r[3]) soilm_w[idx] = soild;
                        solp.zend = zend_r[k]; solp.zenr = zenr_r[k];
                        kstruct kpp = cankCpp(zenr_r[k], x_v[idxl], si);
                        tsdirstruct tspdir = twostreamdirCpp(tir.pait, tir.om, tir.a, tir.gma, tir.J, tir.del, tir.h,
                            gref_m(i,j), kpp.kd, tir.u1, tir.S1, tir.D1, tir.D2);
                        radmodel2 radm = twostreamCpp(pai_v[idxl], clump_v[idxl], gref_m(i,j), svfa_m(i,j), si,
                            tc_v[k], Rsw_v[k], Rdif_v[k], lwdown_v[k], solp, kpp, tspdir, tir);
                        radCsw[hr] = radm.radCsw; radClw[hr] = radm.radClw;
                        Rddown[hr] = radm.Rddown; Rbdown[hr] = radm.Rbdown;
                        radLsw[hr] = radm.radLsw; radLpar[hr] = radm.radLpar;
                        if (out_r[5]) Rdirdown_w[idx] = radm.Rbdown;
                        if (out_r[6]) Rdifdown_w[idx] = radm.Rddown;
                        if (out_r[8]) Rswup_w[idx] = radm.Rdup;
                        double reqhgt2 = (reqhgt < 0.00001) ? 0.00001 : reqhgt;
                        windmodel windm = windCpp(reqhgt2, zref, hgt_v[idxl], pai_v[idxl], u2_v[k], umu_v[k], ws, tiw);
                        uf[hr] = windm.uf; uzday[hr] = windm.uz; gHa[hr] = windm.gHa;
                        if (out_r[4]) uz_w[idx] = windm.uz;
                        soilmodelG0 sG0 = soiltempG0(tc_v[k], es_v[k], ea_v[k], pk_v[k], radm.radGsw, radm.radGlw,
                            tdew_v[k], windm.gHa, soild, mxtc, spa);
                        double Rval = std::abs(sG0.Rnet);
                        if (Rmx < Rval) Rmx = Rval;
                        if (tmx < sG0.Tg) tmx = sG0.Tg;
                        if (tmn > sG0.Tg) tmn = sG0.Tg;
                        surfwet[hr] = sG0.surfwet; radabs[hr] = sG0.radabs;
                    }
                    double dtr = tmx - tmn;
                    for (int hr = 0; hr < 24; ++hr) {
                        int k = dy * 24 + hr + st_r[lyr];
                        int idx = base + cols * rows * k;
                        soilmodel Gvars = soiltemp_hrCpp(tc_v[k], es_v[k], ea_v[k], pk_v[k], radabs[hr], surfwet[hr],
                            tdew_v[k], gHa[hr], soilmday[hr], mxtc, Gp_v[k], dtr, dtrp_v[k], muGp_v[k], kp_v[k], Rmx, sp, spa);
                        Tg[k] = Gvars.Tg; DD[k] = Gvars.DD;
                        if (reqhgt >= 0.0) {
                            radmodel2 rvars;
                            rvars.zend = zend_r[k]; rvars.radCsw = radCsw[hr]; rvars.radClw = radClw[hr];
                            rvars.Rddown = Rddown[hr]; rvars.Rbdown = Rbdown[hr];
                            rvars.radLsw = radLsw[hr]; rvars.radLpar = radLpar[hr];
                            windmodel wvars; wvars.uf = uf[hr]; wvars.uz = uzday[hr]; wvars.gHa = gHa[hr];
                            double reqhgt2 = (reqhgt < 0.00001) ? 0.00001 : reqhgt;
                            abovemodel tv = TVaboveground(reqhgt2, zref, tc_v[k], pk_v[k], ea_v[k], es_v[k], tdew_v[k],
                                Rsw_v[k], Rdif_v[k], lwdown_v[k], soilmday[hr], hgt_v[idxl], pai_v[idxl], paia_v[idxl],
                                x_v[idxl], leafd_v[idxl], leafden_v[idxl], Smin_m(i,j), Smax_m(i,j), Psie_m(i,j),
                                soilb_m(i,j), gsmax_v[idxl], mxtc, stomp, tir, rvars, tiw, wvars, Gvars);
                            if (reqhgt > 0.0) { if (out_r[0]) Tz_w[idx] = tv.Tz; }
                            else              { if (out_r[0]) Tz_w[idx] = Tg[k]; }
                            if (out_r[7]) Rlwdown_w[idx] = tv.lwdn;
                            if (out_r[9]) Rlwup_w[idx] = tv.lwup;
                            if (reqhgt > 0.0) {
                                if (out_r[1]) tleaf_w[idx] = tv.tleaf;
                                if (out_r[2]) relhum_w[idx] = tv.rh;
                            }
                        }
                    }
                }
            }
            if (reqhgt < 0.0 && out_r[0]) {
                double sumD = 0.0;
                for (int k = 0; k < tsteps; ++k) sumD += DD[k];
                double meanD = sumD / static_cast<double>(tsteps);
                std::vector<double> Tzv = Tbelowgroundv(reqhgt, Tg, Tgp2_r, Tbp2_r, meanD, mat, hiy, complete);
                for (int k = 0; k < tsteps; ++k)
                    Tz_w[base + cols * rows * k] = Tzv[k];
            }
        }
    }
};

// [[Rcpp::export]]
List runmicro3Par(DataFrame dfsel, DataFrame obstime, DataFrame climdata, DataFrame pointm,
    List vegp, List soilc, double reqhgt, double zref, double lat, double lon,
    double Sminp, double Smaxp, double tfact, bool complete, double mat,
    std::vector<bool> out, int ncores = 2)
{
    IntegerVector lyr = dfsel["lyr"]; IntegerVector st = dfsel["st"]; IntegerVector ed = dfsel["ed"];
    int nlyrs = lyr.size();
    std::vector<int> ndays_v(nlyrs);
    for (int i = 0; i < nlyrs; ++i) ndays_v[i] = (ed[i] - st[i] + 1) / 24;
    IntegerVector year = obstime["year"]; IntegerVector month = obstime["month"];
    IntegerVector day = obstime["day"];   NumericVector hour = obstime["hour"];
    NumericVector tc = climdata["temp"]; NumericVector es = climdata["es"];
    NumericVector ea = climdata["ea"];   NumericVector tdew = climdata["tdew"];
    NumericVector pk = climdata["pres"]; NumericVector Rsw = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"]; NumericVector lwdown = climdata["lwdown"];
    NumericVector u2 = climdata["windspeed"]; NumericVector wdir = climdata["winddir"];
    NumericVector soilmp = pointm["soilm"]; NumericVector Tgp = pointm["Tg"];
    NumericVector Tbp = pointm["Tbp"];      NumericVector Gp = pointm["G"];
    NumericVector umu = pointm["umu"]; NumericVector kp = pointm["kp"];
    NumericVector muGp = pointm["muGp"]; NumericVector dtrp = pointm["dtrp"];
    NumericVector hgt = vegp["hgt"]; NumericVector pai = vegp["pai"];
    NumericVector x = vegp["x"]; NumericVector gsmax = vegp["gsmax"];
    NumericVector lref = vegp["leafr"]; NumericVector ltra = vegp["leaft"];
    NumericVector clump = vegp["clump"]; NumericVector leafd = vegp["leafd"];
    NumericVector paia = vegp["paia"]; NumericVector leafden = vegp["leafden"];
    NumericMatrix Smin = soilc["Smin"]; NumericMatrix Smax = soilc["Smax"];
    NumericMatrix gref = soilc["gref"]; NumericMatrix soilb = soilc["soilb"];
    NumericMatrix Psie = soilc["Psie"]; NumericMatrix Vq = soilc["Vq"];
    NumericMatrix Vm = soilc["Vm"]; NumericMatrix Mc = soilc["Mc"]; NumericMatrix rho = soilc["rho"];
    NumericMatrix slope = soilc["slope"]; NumericMatrix aspect = soilc["aspect"];
    NumericMatrix twi = soilc["twi"]; NumericMatrix svfa = soilc["svfa"];
    NumericVector wsa = soilc["wsa"]; NumericVector hor = soilc["hor"];
    int rows = Vq.nrow(); int cols = Vq.ncol();
    int tsteps = year.size();
    IntegerVector dim = {rows, cols, tsteps}; int n = rows * cols * tsteps;
    std::vector<double> zend_v(tsteps), zenr_v(tsteps), azid_v(tsteps);
    std::vector<int> sindex_v(tsteps), windex_v(tsteps);
    double mxtc = -273.15;
    for (int k = 0; k < tsteps; ++k) {
        solmodel sp = solpositionCpp(lat, lon, year[k], month[k], day[k], hour[k]);
        zend_v[k] = sp.zend; zenr_v[k] = sp.zenr; azid_v[k] = sp.azid;
        sindex_v[k] = static_cast<int>(std::round(sp.azid / 15.0)) % 24;
        windex_v[k] = static_cast<int>(std::round(wdir[k] / 45)) % 8;
        if (tc[k] > mxtc) mxtc = tc[k];
    }
    int hiy = (year[0] % 4 == 0) ? 366 * 24 : 365 * 24;
    NumericMatrix tadd = soildCppm(twi, Sminp, Smaxp, tfact);
    std::vector<double> Tgp2 = as<std::vector<double>>(Tgp);
    std::vector<double> Tbp2 = as<std::vector<double>>(Tbp);
    std::vector<int> lyr_v(nlyrs), st_v(nlyrs);
    for (int i = 0; i < nlyrs; ++i) { lyr_v[i] = lyr[i]; st_v[i] = st[i]; }
    // Only allocate the output arrays actually requested in `out`. All worker
    // writes are guarded by out_r[k], so unrequested outputs (left empty) are
    // never indexed. Mirrors serial runmicro1-4Cpp; cuts peak output memory by
    // up to 10x when few layers are requested.
    NumericVector Tz, tleaf, relhum, soilmv, uz, Rdirdown, Rdifdown, Rlwdown, Rswup, Rlwup;
    if (out[0]) Tz       = NumericVector(n, NA_REAL);
    if (out[1]) tleaf    = NumericVector(n, NA_REAL);
    if (out[2]) relhum   = NumericVector(n, NA_REAL);
    if (out[3]) soilmv   = NumericVector(n, NA_REAL);
    if (out[4]) uz       = NumericVector(n, NA_REAL);
    if (out[5]) Rdirdown = NumericVector(n, NA_REAL);
    if (out[6]) Rdifdown = NumericVector(n, NA_REAL);
    if (out[7]) Rlwdown  = NumericVector(n, NA_REAL);
    if (out[8]) Rswup    = NumericVector(n, NA_REAL);
    if (out[9]) Rlwup    = NumericVector(n, NA_REAL);
    Micro3Worker worker(
        hgt, pai, x, gsmax, lref, ltra, clump, leafd, paia, leafden,
        Smin, Smax, gref, soilb, Psie, Vq, Vm, Mc, rho, slope, aspect, svfa, tadd,
        wsa, hor, tc, es, ea, tdew, pk, Rsw, Rdif, lwdown, u2,
        soilmp, Gp, umu, kp, muGp, dtrp,
        zend_v, zenr_v, azid_v, sindex_v, windex_v, Tgp2, Tbp2, out,
        lyr_v, st_v, ndays_v, nlyrs,
        reqhgt, zref, lat, mxtc, mat, rows, cols, tsteps, hiy, complete,
        Tz, tleaf, relhum, soilmv, uz, Rdirdown, Rdifdown, Rlwdown, Rswup, Rlwup);
    parallelFor(0, rows * cols, worker, 1, ncores);
    if (out[0]) Tz.attr("dim") = dim;
    if (out[1]) tleaf.attr("dim") = dim;
    if (out[2]) relhum.attr("dim") = dim;
    if (out[3]) soilmv.attr("dim") = dim;
    if (out[4]) uz.attr("dim") = dim;
    if (out[5]) Rdirdown.attr("dim") = dim;
    if (out[6]) Rdifdown.attr("dim") = dim;
    if (out[7]) Rlwdown.attr("dim") = dim;
    if (out[8]) Rswup.attr("dim") = dim;
    if (out[9]) Rlwup.attr("dim") = dim;
    Rcpp::List outp;
    if (out[0]) outp["Tz"] = Tz;
    if (out[1]) outp["tleaf"] = tleaf;
    if (out[2]) outp["relhum"] = relhum;
    if (out[3]) outp["soilm"] = soilmv;
    if (out[4]) outp["windspeed"] = uz;
    if (out[5]) outp["Rdirdown"] = Rdirdown;
    if (out[6]) outp["Rdifdown"] = Rdifdown;
    if (out[7]) outp["Rlwdown"] = Rlwdown;
    if (out[8]) outp["Rswup"] = Rswup;
    if (out[9]) outp["Rlwup"] = Rlwup;
    return outp;
}

// ============================================================
// Micro4Worker: hourly, time-variant veg, array climate
// ============================================================
struct Micro4Worker : public Worker {
    const RVector<double> hgt_v, pai_v, x_v, gsmax_v, lref_v, ltra_v, clump_v, leafd_v, paia_v, leafden_v;
    const RMatrix<double> Smin_m, Smax_m, gref_m, soilb_m, Psie_m, Vq_m, Vm_m, Mc_m, rho_m;
    const RMatrix<double> slope_m, aspect_m, svfa_m, tadd_m, lats_m, lons_m;
    const RVector<double> wsa_v, hor_v;
    const RVector<double> tc_v, es_v, ea_v, tdew_v, pk_v, Rsw_v, Rdif_v, lwdown_v, u2_v;
    const RVector<double> soilmp_v, Tgp_v, Tbp_v, Gp_v, umu_v, kp_v, muGp_v, dtrp_v;
    std::vector<int>    windex_r, year_r, month_r, day_r;
    std::vector<double> hour_r;
    std::vector<bool>   out_r;
    std::vector<int>    lyr_r, st_r, ndays_r;
    int nlyrs;
    double reqhgt, zref, mat;
    int rows, cols, tsteps, hiy;
    bool complete;
    RVector<double> Tz_w, tleaf_w, relhum_w, soilm_w, uz_w;
    RVector<double> Rdirdown_w, Rdifdown_w, Rlwdown_w, Rswup_w, Rlwup_w;

    Micro4Worker(
        NumericVector hgt, NumericVector pai, NumericVector x, NumericVector gsmax,
        NumericVector lref, NumericVector ltra, NumericVector clump, NumericVector leafd,
        NumericVector paia, NumericVector leafden,
        NumericMatrix Smin, NumericMatrix Smax, NumericMatrix gref, NumericMatrix soilb,
        NumericMatrix Psie, NumericMatrix Vq, NumericMatrix Vm, NumericMatrix Mc, NumericMatrix rho,
        NumericMatrix slope, NumericMatrix aspect, NumericMatrix svfa, NumericMatrix tadd,
        NumericMatrix lats, NumericMatrix lons,
        NumericVector wsa, NumericVector hor,
        NumericVector tc, NumericVector es, NumericVector ea, NumericVector tdew,
        NumericVector pk, NumericVector Rsw, NumericVector Rdif, NumericVector lwdown, NumericVector u2,
        NumericVector soilmp, NumericVector Tgp, NumericVector Tbp, NumericVector Gp,
        NumericVector umu, NumericVector kp, NumericVector muGp, NumericVector dtrp,
        std::vector<int> windex, std::vector<int> year, std::vector<int> month,
        std::vector<int> day, std::vector<double> hour,
        std::vector<bool> out,
        std::vector<int> lyr, std::vector<int> st, std::vector<int> ndays, int nlyrs_,
        double reqhgt_, double zref_, double mat_,
        int rows_, int cols_, int tsteps_, int hiy_, bool complete_,
        NumericVector Tz, NumericVector tleaf, NumericVector relhum, NumericVector soilm,
        NumericVector uz, NumericVector Rdirdown, NumericVector Rdifdown,
        NumericVector Rlwdown, NumericVector Rswup, NumericVector Rlwup
    ) :
        hgt_v(hgt), pai_v(pai), x_v(x), gsmax_v(gsmax), lref_v(lref), ltra_v(ltra),
        clump_v(clump), leafd_v(leafd), paia_v(paia), leafden_v(leafden),
        Smin_m(Smin), Smax_m(Smax), gref_m(gref), soilb_m(soilb), Psie_m(Psie),
        Vq_m(Vq), Vm_m(Vm), Mc_m(Mc), rho_m(rho),
        slope_m(slope), aspect_m(aspect), svfa_m(svfa), tadd_m(tadd), lats_m(lats), lons_m(lons),
        wsa_v(wsa), hor_v(hor),
        tc_v(tc), es_v(es), ea_v(ea), tdew_v(tdew), pk_v(pk),
        Rsw_v(Rsw), Rdif_v(Rdif), lwdown_v(lwdown), u2_v(u2),
        soilmp_v(soilmp), Tgp_v(Tgp), Tbp_v(Tbp), Gp_v(Gp),
        umu_v(umu), kp_v(kp), muGp_v(muGp), dtrp_v(dtrp),
        windex_r(std::move(windex)), year_r(std::move(year)), month_r(std::move(month)),
        day_r(std::move(day)), hour_r(std::move(hour)), out_r(std::move(out)),
        lyr_r(std::move(lyr)), st_r(std::move(st)), ndays_r(std::move(ndays)), nlyrs(nlyrs_),
        reqhgt(reqhgt_), zref(zref_), mat(mat_),
        rows(rows_), cols(cols_), tsteps(tsteps_), hiy(hiy_), complete(complete_),
        Tz_w(Tz), tleaf_w(tleaf), relhum_w(relhum), soilm_w(soilm), uz_w(uz),
        Rdirdown_w(Rdirdown), Rdifdown_w(Rdifdown), Rlwdown_w(Rlwdown),
        Rswup_w(Rswup), Rlwup_w(Rlwup)
    {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t cell = begin; cell < end; ++cell) {
            int i = static_cast<int>(cell) % rows;
            int j = static_cast<int>(cell) / rows;
            int base = i + rows * j;
            if (std::isnan(hgt_v[base])) continue;
            soilpstruct spa;
            spa.Smax = Smax_m(i,j); spa.Smin = Smin_m(i,j); spa.soilb = soilb_m(i,j); spa.psi_e = Psie_m(i,j);
            spa.Vq = Vq_m(i,j); spa.Vm = Vm_m(i,j); spa.Mc = Mc_m(i,j); spa.rho = rho_m(i,j);
            soilstruct sp = soilpfun(Vm_m(i,j), Vq_m(i,j), Mc_m(i,j), rho_m(i,j));
            double mxtc = -273.15;
            for (int k = 0; k < tsteps; ++k) {
                double v = tc_v[base + cols * rows * k];
                if (v > mxtc) mxtc = v;
            }
            std::vector<double> Tg(tsteps), DD(tsteps);
            for (int lyr = 0; lyr < nlyrs; ++lyr) {
                int idxl = base + cols * rows * lyr;
                tirstruct tir = twostreamdif(pai_v[idxl], paia_v[idxl], x_v[idxl], lref_v[idxl], ltra_v[idxl], clump_v[idxl], gref_m(i,j));
                stompstruct stomp = stomparamsCpp(hgt_v[idxl], lats_m(i,j), x_v[idxl]);
                tiwstruct tiw = windtiCpp(hgt_v[idxl], pai_v[idxl]);
                for (int dy = 0; dy < ndays_r[lyr]; ++dy) {
                    double Rmx = -999.9, tmx = -999.0, tmn = 999.0;
                    std::vector<double> surfwet(24), radabs(24), soilmday(24);
                    std::vector<double> radCsw(24), radClw(24), Rddown(24), Rbdown(24), radLsw(24), radLpar(24);
                    std::vector<double> uf(24), uzday(24), gHa(24), zendday(24);
                    for (int hr = 0; hr < 24; ++hr) {
                        int k = dy * 24 + hr + st_r[lyr];
                        int idx = base + cols * rows * k;
                        solmodel solp = solpositionCpp(lats_m(i,j), lons_m(i,j), year_r[k], month_r[k], day_r[k], hour_r[k]);
                        zendday[hr] = solp.zend;
                        double si = solarindexCpp(slope_m(i,j), aspect_m(i,j), solp.zend, solp.azid);
                        if (si < 0.0) si = 0.0;
                        int sindx = static_cast<int>(std::round(solp.azid / 15)) % 24;
                        double ha = hor_v[sindx * rows * cols + j * rows + i];
                        double sa = (pi / 2.0) - solp.zenr;
                        if (ha > std::tan(sa)) si = 0.0;
                        double ws = wsa_v[windex_r[k] * rows * cols + j * rows + i];
                        double soild = soildCpp(soilmp_v[idx], Smin_m(i,j), Smax_m(i,j), tadd_m(i,j));
                        soilmday[hr] = soild;
                        if (out_r[3]) soilm_w[idx] = soild;
                        kstruct kpp = cankCpp(solp.zenr, x_v[idxl], si);
                        tsdirstruct tspdir = twostreamdirCpp(tir.pait, tir.om, tir.a, tir.gma, tir.J, tir.del, tir.h,
                            gref_m(i,j), kpp.kd, tir.u1, tir.S1, tir.D1, tir.D2);
                        radmodel2 radm = twostreamCpp(pai_v[idxl], clump_v[idxl], gref_m(i,j), svfa_m(i,j), si,
                            tc_v[idx], Rsw_v[idx], Rdif_v[idx], lwdown_v[idx], solp, kpp, tspdir, tir);
                        radCsw[hr] = radm.radCsw; radClw[hr] = radm.radClw;
                        Rddown[hr] = radm.Rddown; Rbdown[hr] = radm.Rbdown;
                        radLsw[hr] = radm.radLsw; radLpar[hr] = radm.radLpar;
                        if (out_r[5]) Rdirdown_w[idx] = radm.Rbdown;
                        if (out_r[6]) Rdifdown_w[idx] = radm.Rddown;
                        if (out_r[8]) Rswup_w[idx] = radm.Rdup;
                        double reqhgt2 = (reqhgt < 0.00001) ? 0.00001 : reqhgt;
                        windmodel windm = windCpp(reqhgt2, zref, hgt_v[idxl], pai_v[idxl], u2_v[idx], umu_v[idx], ws, tiw);
                        uf[hr] = windm.uf; uzday[hr] = windm.uz; gHa[hr] = windm.gHa;
                        if (out_r[4]) uz_w[idx] = windm.uz;
                        soilmodelG0 sG0 = soiltempG0(tc_v[idx], es_v[idx], ea_v[idx], pk_v[idx], radm.radGsw, radm.radGlw,
                            tdew_v[idx], windm.gHa, soild, mxtc, spa);
                        double Rval = std::abs(sG0.Rnet);
                        if (Rmx < Rval) Rmx = Rval;
                        if (tmx < sG0.Tg) tmx = sG0.Tg;
                        if (tmn > sG0.Tg) tmn = sG0.Tg;
                        surfwet[hr] = sG0.surfwet; radabs[hr] = sG0.radabs;
                    }
                    double dtr = tmx - tmn;
                    for (int hr = 0; hr < 24; ++hr) {
                        int k = dy * 24 + hr + st_r[lyr];
                        int idx = base + cols * rows * k;
                        soilmodel Gvars = soiltemp_hrCpp(tc_v[idx], es_v[idx], ea_v[idx], pk_v[idx], radabs[hr], surfwet[hr],
                            tdew_v[idx], gHa[hr], soilmday[hr], mxtc, Gp_v[idx], dtr, dtrp_v[idx], muGp_v[idx], kp_v[idx], Rmx, sp, spa);
                        Tg[k] = Gvars.Tg; DD[k] = Gvars.DD;
                        if (reqhgt >= 0.0) {
                            radmodel2 rvars;
                            rvars.zend = zendday[hr]; rvars.radCsw = radCsw[hr]; rvars.radClw = radClw[hr];
                            rvars.Rddown = Rddown[hr]; rvars.Rbdown = Rbdown[hr];
                            rvars.radLsw = radLsw[hr]; rvars.radLpar = radLpar[hr];
                            windmodel wvars; wvars.uf = uf[hr]; wvars.uz = uzday[hr]; wvars.gHa = gHa[hr];
                            double reqhgt2 = (reqhgt < 0.00001) ? 0.00001 : reqhgt;
                            abovemodel tv = TVaboveground(reqhgt2, zref, tc_v[idx], pk_v[idx], ea_v[idx], es_v[idx], tdew_v[idx],
                                Rsw_v[idx], Rdif_v[idx], lwdown_v[idx], soilmday[hr], hgt_v[idxl], pai_v[idxl], paia_v[idxl],
                                x_v[idxl], leafd_v[idxl], leafden_v[idxl], Smin_m(i,j), Smax_m(i,j), Psie_m(i,j),
                                soilb_m(i,j), gsmax_v[idxl], mxtc, stomp, tir, rvars, tiw, wvars, Gvars);
                            if (reqhgt > 0.0) { if (out_r[0]) Tz_w[idx] = tv.Tz; }
                            else              { if (out_r[0]) Tz_w[idx] = Tg[k]; }
                            if (out_r[7]) Rlwdown_w[idx] = tv.lwdn;
                            if (out_r[9]) Rlwup_w[idx] = tv.lwup;
                            if (reqhgt > 0.0) {
                                if (out_r[1]) tleaf_w[idx] = tv.tleaf;
                                if (out_r[2]) relhum_w[idx] = tv.rh;
                            }
                        }
                    }
                }
            }
            if (reqhgt < 0.0 && out_r[0]) {
                std::vector<double> Tgpv(tsteps), Tbpv(tsteps);
                double sumD = 0.0;
                for (int k = 0; k < tsteps; ++k) {
                    int idx = base + cols * rows * k;
                    Tgpv[k] = Tgp_v[idx]; Tbpv[k] = Tbp_v[idx];
                    sumD += DD[k];
                }
                double meanD = sumD / static_cast<double>(tsteps);
                std::vector<double> Tzv = Tbelowgroundv(reqhgt, Tg, Tgpv, Tbpv, meanD, mat, hiy, complete);
                for (int k = 0; k < tsteps; ++k)
                    Tz_w[base + cols * rows * k] = Tzv[k];
            }
        }
    }
};

// [[Rcpp::export]]
List runmicro4Par(DataFrame dfsel, DataFrame obstime, List climdata, List pointm, List vegp,
    List soilc, double reqhgt, double zref, NumericMatrix lats, NumericMatrix lons,
    double Sminp, double Smaxp, double tfact, bool complete, double mat,
    std::vector<bool> out, int ncores = 2)
{
    IntegerVector lyr = dfsel["lyr"]; IntegerVector st = dfsel["st"]; IntegerVector ed = dfsel["ed"];
    int nlyrs = lyr.size();
    std::vector<int> ndays_v(nlyrs);
    for (int i = 0; i < nlyrs; ++i) ndays_v[i] = (ed[i] - st[i] + 1) / 24;
    IntegerVector year = obstime["year"]; IntegerVector month = obstime["month"];
    IntegerVector day = obstime["day"];   NumericVector hour = obstime["hour"];
    NumericVector tc = climdata["tc"]; NumericVector es = climdata["es"];
    NumericVector ea = climdata["ea"]; NumericVector tdew = climdata["tdew"];
    NumericVector pk = climdata["pk"]; NumericVector Rsw = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"]; NumericVector lwdown = climdata["lwdown"];
    NumericVector u2 = climdata["windspeed"]; NumericVector wdir = climdata["winddir"];
    NumericVector soilmp = pointm["soilm"]; NumericVector Tgp = pointm["Tg"];
    NumericVector Tbp = pointm["Tbp"]; NumericVector Gp = pointm["Gp"];
    NumericVector umu = pointm["umu"]; NumericVector kp = pointm["kp"];
    NumericVector muGp = pointm["muGp"]; NumericVector dtrp = pointm["dtrp"];
    NumericVector hgt = vegp["hgt"]; NumericVector pai = vegp["pai"];
    NumericVector x = vegp["x"]; NumericVector gsmax = vegp["gsmax"];
    NumericVector lref = vegp["leafr"]; NumericVector ltra = vegp["leaft"];
    NumericVector clump = vegp["clump"]; NumericVector leafd = vegp["leafd"];
    NumericVector paia = vegp["paia"]; NumericVector leafden = vegp["leafden"];
    NumericMatrix Smin = soilc["Smin"]; NumericMatrix Smax = soilc["Smax"];
    NumericMatrix gref = soilc["gref"]; NumericMatrix soilb = soilc["soilb"];
    NumericMatrix Psie = soilc["Psie"]; NumericMatrix Vq = soilc["Vq"];
    NumericMatrix Vm = soilc["Vm"]; NumericMatrix Mc = soilc["Mc"]; NumericMatrix rho = soilc["rho"];
    NumericMatrix slope = soilc["slope"]; NumericMatrix aspect = soilc["aspect"];
    NumericMatrix twi = soilc["twi"]; NumericMatrix svfa = soilc["svfa"];
    NumericVector wsa = soilc["wsa"]; NumericVector hor = soilc["hor"];
    int rows = Smin.nrow(); int cols = Smin.ncol();
    int tsteps = year.size();
    IntegerVector dim = {rows, cols, tsteps}; int n = rows * cols * tsteps;
    std::vector<int> windex_v(tsteps);
    for (int k = 0; k < tsteps; ++k)
        windex_v[k] = static_cast<int>(std::round(wdir[k] / 45)) % 8;
    int hiy = (year[0] % 4 == 0) ? 366 * 24 : 365 * 24;
    NumericMatrix tadd = soildCppm(twi, Sminp, Smaxp, tfact);
    std::vector<int> yv(tsteps), mv(tsteps), dv(tsteps);
    std::vector<double> hv(tsteps);
    for (int k = 0; k < tsteps; ++k) { yv[k]=year[k]; mv[k]=month[k]; dv[k]=day[k]; hv[k]=hour[k]; }
    std::vector<int> lyr_v(nlyrs), st_v(nlyrs);
    for (int i = 0; i < nlyrs; ++i) { lyr_v[i] = lyr[i]; st_v[i] = st[i]; }
    // Only allocate the output arrays actually requested in `out`. All worker
    // writes are guarded by out_r[k], so unrequested outputs (left empty) are
    // never indexed. Mirrors serial runmicro1-4Cpp; cuts peak output memory by
    // up to 10x when few layers are requested.
    NumericVector Tz, tleaf, relhum, soilmv, uz, Rdirdown, Rdifdown, Rlwdown, Rswup, Rlwup;
    if (out[0]) Tz       = NumericVector(n, NA_REAL);
    if (out[1]) tleaf    = NumericVector(n, NA_REAL);
    if (out[2]) relhum   = NumericVector(n, NA_REAL);
    if (out[3]) soilmv   = NumericVector(n, NA_REAL);
    if (out[4]) uz       = NumericVector(n, NA_REAL);
    if (out[5]) Rdirdown = NumericVector(n, NA_REAL);
    if (out[6]) Rdifdown = NumericVector(n, NA_REAL);
    if (out[7]) Rlwdown  = NumericVector(n, NA_REAL);
    if (out[8]) Rswup    = NumericVector(n, NA_REAL);
    if (out[9]) Rlwup    = NumericVector(n, NA_REAL);
    Micro4Worker worker(
        hgt, pai, x, gsmax, lref, ltra, clump, leafd, paia, leafden,
        Smin, Smax, gref, soilb, Psie, Vq, Vm, Mc, rho, slope, aspect, svfa, tadd, lats, lons,
        wsa, hor, tc, es, ea, tdew, pk, Rsw, Rdif, lwdown, u2,
        soilmp, Tgp, Tbp, Gp, umu, kp, muGp, dtrp,
        windex_v, yv, mv, dv, hv, out,
        lyr_v, st_v, ndays_v, nlyrs,
        reqhgt, zref, mat, rows, cols, tsteps, hiy, complete,
        Tz, tleaf, relhum, soilmv, uz, Rdirdown, Rdifdown, Rlwdown, Rswup, Rlwup);
    parallelFor(0, rows * cols, worker, 1, ncores);
    if (out[0]) Tz.attr("dim") = dim;
    if (out[1]) tleaf.attr("dim") = dim;
    if (out[2]) relhum.attr("dim") = dim;
    if (out[3]) soilmv.attr("dim") = dim;
    if (out[4]) uz.attr("dim") = dim;
    if (out[5]) Rdirdown.attr("dim") = dim;
    if (out[6]) Rdifdown.attr("dim") = dim;
    if (out[7]) Rlwdown.attr("dim") = dim;
    if (out[8]) Rswup.attr("dim") = dim;
    if (out[9]) Rlwup.attr("dim") = dim;
    Rcpp::List outp;
    if (out[0]) outp["Tz"] = Tz;
    if (out[1]) outp["tleaf"] = tleaf;
    if (out[2]) outp["relhum"] = relhum;
    if (out[3]) outp["soilm"] = soilmv;
    if (out[4]) outp["windspeed"] = uz;
    if (out[5]) outp["Rdirdown"] = Rdirdown;
    if (out[6]) outp["Rdifdown"] = Rdifdown;
    if (out[7]) outp["Rlwdown"] = Rlwdown;
    if (out[8]) outp["Rswup"] = Rswup;
    if (out[9]) outp["Rlwup"] = Rlwup;
    return outp;
}

// ============================================================
// Snow1Worker: snow model, data.frame climate
// ============================================================
struct Snow1Worker : public Worker {
    const RMatrix<double> pai_m, hgt_m, ltra_m, clump_m;
    const RMatrix<double> slope_m, aspect_m, skyview_m;
    const RMatrix<double> isnowdc_m, isnowdg_m;
    const RMatrix<int>    isnowac_m, isnowag_m;
    const RVector<double> wsa_v, hor_v;
    const RVector<double> tc_v, pk_v, Rsw_v, Rdif_v, Rlw_v, u2_v, prec_v;
    const RVector<double> ea_v, te_v, Gp_v, Tcp_v, umu_v;
    const RVector<double> Rmx_v, Rmn_v, Rswmx_v, Rlwmx_v, Rswmn_v, Rlwmn_v, Gmx_v, salb_v;
    std::vector<double> zend_r, azid_r;
    std::vector<int>    sindex_r, windex_r;
    std::vector<int>    year_r, month_r, day_r;
    std::vector<double> hour_r;
    std::vector<double> sdp_r;
    double lat, lon, zref;
    int rows, cols, tsteps;
    RVector<double> Tc_w, Tg_w, sdepc_w, sdepg_w, sdenc_w;
    RMatrix<double> agec_w, ageg_w, meltc_w, meltg_w;

    Snow1Worker(
        NumericMatrix pai, NumericMatrix hgt, NumericMatrix ltra, NumericMatrix clump,
        NumericMatrix slope, NumericMatrix aspect, NumericMatrix skyview,
        NumericMatrix isnowdc, NumericMatrix isnowdg,
        IntegerMatrix isnowac, IntegerMatrix isnowag,
        NumericVector wsa, NumericVector hor,
        NumericVector tc, NumericVector pk, NumericVector Rsw, NumericVector Rdif,
        NumericVector Rlw, NumericVector u2, NumericVector prec,
        NumericVector ea, NumericVector te, NumericVector Gp, NumericVector Tcp, NumericVector umu,
        NumericVector Rmx, NumericVector Rmn, NumericVector Rswmx, NumericVector Rlwmx,
        NumericVector Rswmn, NumericVector Rlwmn, NumericVector Gmx, NumericVector salb,
        std::vector<double> zend, std::vector<double> azid,
        std::vector<int> sindex, std::vector<int> windex,
        std::vector<int> year, std::vector<int> month, std::vector<int> day, std::vector<double> hour,
        std::vector<double> sdp,
        double lat_, double lon_, double zref_,
        int rows_, int cols_, int tsteps_,
        NumericVector Tc, NumericVector Tg, NumericVector sdepc, NumericVector sdepg, NumericVector sdenc,
        NumericMatrix agec, NumericMatrix ageg, NumericMatrix meltc, NumericMatrix meltg
    ) :
        pai_m(pai), hgt_m(hgt), ltra_m(ltra), clump_m(clump),
        slope_m(slope), aspect_m(aspect), skyview_m(skyview),
        isnowdc_m(isnowdc), isnowdg_m(isnowdg), isnowac_m(isnowac), isnowag_m(isnowag),
        wsa_v(wsa), hor_v(hor),
        tc_v(tc), pk_v(pk), Rsw_v(Rsw), Rdif_v(Rdif), Rlw_v(Rlw), u2_v(u2), prec_v(prec),
        ea_v(ea), te_v(te), Gp_v(Gp), Tcp_v(Tcp), umu_v(umu),
        Rmx_v(Rmx), Rmn_v(Rmn), Rswmx_v(Rswmx), Rlwmx_v(Rlwmx),
        Rswmn_v(Rswmn), Rlwmn_v(Rlwmn), Gmx_v(Gmx), salb_v(salb),
        zend_r(std::move(zend)), azid_r(std::move(azid)),
        sindex_r(std::move(sindex)), windex_r(std::move(windex)),
        year_r(std::move(year)), month_r(std::move(month)), day_r(std::move(day)), hour_r(std::move(hour)),
        sdp_r(std::move(sdp)),
        lat(lat_), lon(lon_), zref(zref_), rows(rows_), cols(cols_), tsteps(tsteps_),
        Tc_w(Tc), Tg_w(Tg), sdepc_w(sdepc), sdepg_w(sdepg), sdenc_w(sdenc),
        agec_w(agec), ageg_w(ageg), meltc_w(meltc), meltg_w(meltg)
    {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t cell = begin; cell < end; ++cell) {
            int i = static_cast<int>(cell) % rows;
            int j = static_cast<int>(cell) / rows;
            if (std::isnan(hgt_m(i, j))) continue;
            obspoint obstimeo{}; climpoint climo{}; vegpoint vegpo{};
            otherpoint othero{}; snowpoint snowo{};
            vegpo.pai = pai_m(i,j); vegpo.hgt = hgt_m(i,j); vegpo.clump = clump_m(i,j); vegpo.ltra = ltra_m(i,j);
            othero.slope = slope_m(i,j); othero.aspect = aspect_m(i,j);
            othero.lat = lat; othero.lon = lon; othero.zref = zref; othero.psim = 0.0; othero.psih = 0.0;
            int snowagec = isnowac_m(i, j);
            int snowageg = isnowag_m(i, j);
            double sdencp = ((sdp_r[0] - sdp_r[1]) * (1 - std::exp(-sdp_r[2] * isnowdc_m(i,j) / 100.0 -
                sdp_r[3] * snowagec / 24.0)) + sdp_r[1]) * 1000.0;
            double sdengp = ((sdp_r[0] - sdp_r[1]) * (1 - std::exp(-sdp_r[2] * isnowdg_m(i,j) * 0.5 / 100.0 -
                sdp_r[3] * snowageg / 24.0)) + sdp_r[1]) * 1000.0;
            double sdepcp = isnowdc_m(i,j);
            double sdepgp = isnowdg_m(i,j);
            meltc_w(i, j) = 0.0;
            meltg_w(i, j) = 0.0;
            for (int k = 0; k < tsteps; ++k) {
                int idx = i + rows * j + cols * rows * k;
                int snowtest = 0;
                if (sdepcp > 0.0) snowtest = 1;
                if (tc_v[k] < 2.0 && prec_v[k] > 0.0) snowtest = 1;
                if (snowtest > 0) {
                    double paip = pai_m(i,j);
                    if (hgt_m(i,j) > sdepgp) paip = paip * (hgt_m(i,j) - sdepgp) / hgt_m(i,j);
                    double dtR = Rmx_v[k] - Rmn_v[k];
                    double trS = skyview_m(i,j) * std::exp(-paip);
                    double Rem = 0.97 * sb * radem(tc_v[k]);
                    double dmxS = trS * Rswmx_v[k] + trS * Rlwmx_v[k] + (1 - trS) * Rem - Rem;
                    double dmnS = trS * Rswmn_v[k] + trS * Rlwmn_v[k] + (1 - trS) * Rem - Rem;
                    double Gmu = (dmxS - dmnS) / dtR;
                    double G = Gp_v[k] * Gmu;
                    if (G > Gmx_v[k]) G = Gmx_v[k];
                    if (G < -Gmx_v[k]) G = -Gmx_v[k];
                    double ha = hor_v[sindex_r[k] * rows * cols + j * rows + i];
                    double sa = 90.0 - zend_r[k];
                    double smu = 1.0;
                    if (ha > std::tan(sa * torad)) smu = 0.0;
                    double ws = wsa_v[windex_r[k] * rows * cols + j * rows + i];
                    double u2p = umu_v[k] * ws * u2_v[k];
                    double Rdifp = Rdif_v[k] * skyview_m(i,j);
                    double Rdirp = (Rsw_v[k] - Rdif_v[k]) * smu;
                    double Rswp = Rdirp + Rdifp;
                    double Rlwp = Rlw_v[k] * skyview_m(i,j);
                    obstimeo.year = year_r[k]; obstimeo.month = month_r[k];
                    obstimeo.day = day_r[k];   obstimeo.hour = hour_r[k];
                    climo.tc = tc_v[k]; climo.ea = ea_v[k]; climo.pk = pk_v[k];
                    climo.u2 = u2p; climo.prec = prec_v[k];
                    climo.Rsw = Rswp; climo.Rdif = Rdifp; climo.Rlw = Rlwp;
                    climo.Tci = Tcp_v[k]; climo.te = te_v[k];
                    othero.G = G;
                    snowo.alb = salb_v[k]; snowo.sdenc = sdencp; snowo.sdeng = sdengp;
                    snowo.sdepc = sdepcp; snowo.sdepg = sdepgp;
                    snowo.snowagec = static_cast<double>(snowagec);
                    snowo.snowageg = static_cast<double>(snowageg);
                    snowmodpoint smod = snowoneB(obstimeo, climo, vegpo, snowo, othero, sdp_r, 1.0);
                    Tc_w[idx] = smod.Tc; Tg_w[idx] = smod.Tg;
                    sdepc_w[idx] = smod.sdepc; sdepg_w[idx] = smod.sdepg; sdenc_w[idx] = smod.sdenc;
                    sdencp = smod.sdenc; sdengp = smod.sdeng;
                    sdepcp = smod.sdepc; sdepgp = smod.sdepg;
                    snowagec = smod.snowagec; snowageg = smod.snowageg;
                    double melc = smod.mSc + smod.mMc + smod.mRc;
                    double melg = smod.mSg + smod.mMg + smod.mRg;
                    meltc_w(i,j) = meltc_w(i,j) + (melc * 1000.0) / smod.sdenc;
                    meltg_w(i,j) = meltg_w(i,j) + (melg * 1000.0) / smod.sdeng;
                } else {
                    Tc_w[idx] = 0.0; Tg_w[idx] = 0.0;
                    sdepc_w[idx] = 0.0; sdepg_w[idx] = 0.0;
                    sdenc_w[idx] = sdp_r[1] * 1000.0;
                }
            }
            agec_w(i,j) = static_cast<double>(snowagec);
            ageg_w(i,j) = static_cast<double>(snowageg);
        }
    }
};

// [[Rcpp::export]]
List gridmodelsnow1Par(DataFrame obstime, DataFrame climdata, DataFrame pointm, List vegp,
    List other, std::string snowenv, int ncores = 2)
{
    IntegerVector year = obstime["year"]; IntegerVector month = obstime["month"];
    IntegerVector day = obstime["day"];   NumericVector hour = obstime["hour"];
    NumericVector Gp = pointm["Gp"]; NumericVector Tcp = pointm["Tc"];
    NumericVector RswabsG = pointm["RswabsG"]; NumericVector RlwabsG = pointm["RlwabsG"];
    NumericVector umu = pointm["umu"]; NumericVector tr = pointm["tr"];
    NumericMatrix pai = vegp["pai"]; NumericMatrix hgt = vegp["hgt"];
    NumericMatrix ltra = vegp["leaft"]; NumericMatrix clump = vegp["clump"];
    NumericMatrix slope = other["slope"]; NumericMatrix aspect = other["aspect"];
    NumericMatrix skyview = other["skyview"];
    NumericVector wsa = other["wsa"]; NumericVector hor = other["hor"];
    double lat = other["lat"]; double lon = other["lon"]; double zref = other["zref"];
    NumericMatrix isnowdc = other["isnowdc"]; NumericMatrix isnowdg = other["isnowdg"];
    IntegerMatrix isnowac = other["isnowac"]; IntegerMatrix isnowag = other["isnowag"];
    NumericVector tc = climdata["temp"]; NumericVector rh = climdata["relhum"];
    NumericVector pk = climdata["pres"]; NumericVector Rsw = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"]; NumericVector Rlw = climdata["lwdown"];
    NumericVector u2 = climdata["windspeed"]; NumericVector wdir = climdata["winddir"];
    NumericVector prec = climdata["precip"];
    NumericVector salb = snowalbCpp(prec);
    int rows = pai.nrow(); int cols = pai.ncol();
    int tsteps = tc.size(); int ndays = tr.size() / 24;
    IntegerVector dim = {rows, cols, tsteps};
    NumericVector ea(tsteps), te(tsteps), Rnet(tsteps);
    for (int k = 0; k < tsteps; ++k) {
        ea[k] = satvapCpp(tc[k]) * rh[k] / 100.0;
        te[k] = (Tcp[k] + tc[k]) / 2.0;
        double Rem = 0.97 * sb * radem(tc[k]);
        Rnet[k] = RswabsG[k] + RlwabsG[k] - Rem;
    }
    NumericVector Rmxd(ndays, -1352.0), Rmnd(ndays, 1352.0);
    NumericVector Rswmnd(ndays), Rlwmnd(ndays), Rswmxd(ndays), Rlwmxd(ndays);
    for (int d = 0; d < ndays; ++d) {
        for (int h = 0; h < 24; ++h) {
            int k = d * 24 + h;
            if (Rmxd[d] < Rnet[k]) { Rmxd[d] = Rnet[k]; Rswmxd[d] = Rsw[k]; Rlwmxd[d] = Rlw[k]; }
            if (Rmnd[d] > Rnet[k]) { Rmnd[d] = Rnet[k]; Rswmnd[d] = Rsw[k]; Rlwmnd[d] = Rlw[k]; }
        }
    }
    NumericVector Rmx(tsteps), Rmn(tsteps), Rswmn(tsteps), Rlwmn(tsteps), Rswmx(tsteps), Rlwmx(tsteps);
    NumericVector Gmxd(ndays, 0.0), Gmx(tsteps);
    for (int d = 0; d < ndays; ++d) {
        for (int h = 0; h < 24; ++h) {
            int k = d * 24 + h;
            Rmx[k] = Rmxd[d]; Rmn[k] = Rmnd[d];
            Rswmn[k] = Rswmnd[d]; Rlwmn[k] = Rlwmnd[d];
            Rswmx[k] = Rswmxd[d]; Rlwmx[k] = Rlwmxd[d];
            if (std::abs(Rnet[k]) > Gmxd[d]) Gmxd[d] = std::abs(Rnet[k]);
        }
        for (int h = 0; h < 24; ++h) Gmx[d * 24 + h] = Gmxd[d];
    }
    std::vector<double> zend_v(tsteps), azid_v(tsteps);
    std::vector<int> sindex_v(tsteps), windex_v(tsteps);
    for (int k = 0; k < tsteps; ++k) {
        solmodel sp = solpositionCpp(lat, lon, year[k], month[k], day[k], hour[k]);
        zend_v[k] = sp.zend; azid_v[k] = sp.azid;
        sindex_v[k] = static_cast<int>(std::round(sp.azid / 15)) % 24;
        windex_v[k] = static_cast<int>(std::round(wdir[k] / 45)) % 8;
    }
    std::vector<double> sdp = snowdenp(snowenv);
    std::vector<int> yv(tsteps), mv(tsteps), dv(tsteps); std::vector<double> hv(tsteps);
    for (int k = 0; k < tsteps; ++k) { yv[k]=year[k]; mv[k]=month[k]; dv[k]=day[k]; hv[k]=hour[k]; }
    int n3 = rows * cols * tsteps;
    NumericVector Tc(n3, NA_REAL), Tg(n3, NA_REAL), sdepc(n3, NA_REAL), sdepg(n3, NA_REAL), sdenc(n3, NA_REAL);
    NumericMatrix agec = bioclimfill(rows, cols), ageg = bioclimfill(rows, cols);
    NumericMatrix meltc = bioclimfill(rows, cols), meltg = bioclimfill(rows, cols);
    Tc.attr("dim") = dim; Tg.attr("dim") = dim;
    sdepc.attr("dim") = dim; sdepg.attr("dim") = dim; sdenc.attr("dim") = dim;
    Snow1Worker worker(
        pai, hgt, ltra, clump, slope, aspect, skyview, isnowdc, isnowdg, isnowac, isnowag,
        wsa, hor, tc, pk, Rsw, Rdif, Rlw, u2, prec,
        ea, te, Gp, Tcp, umu,
        Rmx, Rmn, Rswmx, Rlwmx, Rswmn, Rlwmn, Gmx, salb,
        zend_v, azid_v, sindex_v, windex_v, yv, mv, dv, hv, sdp,
        lat, lon, zref, rows, cols, tsteps,
        Tc, Tg, sdepc, sdepg, sdenc, agec, ageg, meltc, meltg);
    parallelFor(0, rows * cols, worker, 1, ncores);
    List out;
    out["Tc"] = Tc; out["Tg"] = Tg; out["sdepc"] = sdepc; out["sdepg"] = sdepg;
    out["sden"] = sdenc; out["agec"] = agec; out["ageg"] = ageg;
    out["meltc"] = meltc; out["meltg"] = meltg;
    return out;
}

// ============================================================
// Snow2Worker: snow model, array climate
// ============================================================
struct Snow2Worker : public Worker {
    const RMatrix<double> pai_m, hgt_m, ltra_m, clump_m;
    const RMatrix<double> slope_m, aspect_m, skyview_m, lats_m, lons_m;
    const RMatrix<double> isnowdc_m, isnowdg_m;
    const RMatrix<int>    isnowac_m, isnowag_m;
    const RVector<double> wsa_v, hor_v;
    const RVector<double> tc_v, rh_v, pk_v, Rsw_v, Rdif_v, Rlw_v, u2_v, prec_v;
    const RVector<double> Gp_v, Tcp_v, RswabsG_v, RlwabsG_v, umu_v;
    std::vector<int>    windex_r, year_r, month_r, day_r;
    std::vector<double> hour_r, sdp_r;
    double zref;
    int rows, cols, tsteps;
    RVector<double> Tc_w, Tg_w, sdepc_w, sdepg_w, sdenc_w;
    RMatrix<double> agec_w, ageg_w, meltc_w, meltg_w;

    Snow2Worker(
        NumericMatrix pai, NumericMatrix hgt, NumericMatrix ltra, NumericMatrix clump,
        NumericMatrix slope, NumericMatrix aspect, NumericMatrix skyview,
        NumericMatrix lats, NumericMatrix lons,
        NumericMatrix isnowdc, NumericMatrix isnowdg,
        IntegerMatrix isnowac, IntegerMatrix isnowag,
        NumericVector wsa, NumericVector hor,
        NumericVector tc, NumericVector rh, NumericVector pk, NumericVector Rsw,
        NumericVector Rdif, NumericVector Rlw, NumericVector u2, NumericVector prec,
        NumericVector Gp, NumericVector Tcp, NumericVector RswabsG, NumericVector RlwabsG, NumericVector umu,
        std::vector<int> windex, std::vector<int> year, std::vector<int> month,
        std::vector<int> day, std::vector<double> hour, std::vector<double> sdp,
        double zref_, int rows_, int cols_, int tsteps_,
        NumericVector Tc, NumericVector Tg, NumericVector sdepc, NumericVector sdepg, NumericVector sdenc,
        NumericMatrix agec, NumericMatrix ageg, NumericMatrix meltc, NumericMatrix meltg
    ) :
        pai_m(pai), hgt_m(hgt), ltra_m(ltra), clump_m(clump),
        slope_m(slope), aspect_m(aspect), skyview_m(skyview), lats_m(lats), lons_m(lons),
        isnowdc_m(isnowdc), isnowdg_m(isnowdg), isnowac_m(isnowac), isnowag_m(isnowag),
        wsa_v(wsa), hor_v(hor),
        tc_v(tc), rh_v(rh), pk_v(pk), Rsw_v(Rsw), Rdif_v(Rdif), Rlw_v(Rlw), u2_v(u2), prec_v(prec),
        Gp_v(Gp), Tcp_v(Tcp), RswabsG_v(RswabsG), RlwabsG_v(RlwabsG), umu_v(umu),
        windex_r(std::move(windex)), year_r(std::move(year)), month_r(std::move(month)),
        day_r(std::move(day)), hour_r(std::move(hour)), sdp_r(std::move(sdp)),
        zref(zref_), rows(rows_), cols(cols_), tsteps(tsteps_),
        Tc_w(Tc), Tg_w(Tg), sdepc_w(sdepc), sdepg_w(sdepg), sdenc_w(sdenc),
        agec_w(agec), ageg_w(ageg), meltc_w(meltc), meltg_w(meltg)
    {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t cell = begin; cell < end; ++cell) {
            int i = static_cast<int>(cell) % rows;
            int j = static_cast<int>(cell) / rows;
            if (std::isnan(hgt_m(i, j))) continue;
            int ndays = tsteps / 24;
            std::vector<double> Rnet(tsteps);
            for (int k = 0; k < tsteps; ++k) {
                int idx = i + rows * j + cols * rows * k;
                double Rem = 0.97 * sb * radem(tc_v[idx]);
                Rnet[k] = RswabsG_v[idx] + RlwabsG_v[idx] - Rem;
            }
            // Inline snowalbCpp (thread-safe: std:: types only, matches serial formula exactly)
            std::vector<int> hs(tsteps, 0);
            for (int k = 1; k < tsteps; ++k) {
                int idx = i + rows * j + cols * rows * k;
                hs[k] = (prec_v[idx] > 0.0) ? 0 : hs[k-1] + 1;
            }
            std::vector<double> salb(tsteps);
            for (int k = 0; k < tsteps; ++k) {
                // hs[k] / 24 is integer division — matches snowalbCpp exactly
                salb[k] = (-9.8740 * std::log(hs[k] / 24) + 78.3434) / 100.0;
                if (salb[k] > 0.95) salb[k] = 0.95;
                if (salb[k] < 0.1) salb[k] = 0.1;
            }
            std::vector<double> Rmxd(ndays, -1352.0), Rmnd(ndays, 1352.0);
            std::vector<double> Rswmnd(ndays), Rlwmnd(ndays), Rswmxd(ndays), Rlwmxd(ndays);
            for (int d = 0; d < ndays; ++d) {
                for (int h = 0; h < 24; ++h) {
                    int k = d * 24 + h;
                    int idx = i + rows * j + cols * rows * k;
                    if (Rmxd[d] < Rnet[k]) { Rmxd[d] = Rnet[k]; Rswmxd[d] = Rsw_v[idx]; Rlwmxd[d] = Rlw_v[idx]; }
                    if (Rmnd[d] > Rnet[k]) { Rmnd[d] = Rnet[k]; Rswmnd[d] = Rsw_v[idx]; Rlwmnd[d] = Rlw_v[idx]; }
                }
            }
            std::vector<double> Rmx(tsteps), Rmn(tsteps), Rswmn(tsteps), Rlwmn(tsteps), Rswmx(tsteps), Rlwmx(tsteps);
            std::vector<double> Gmxd(ndays, 0.0), Gmx(tsteps);
            for (int d = 0; d < ndays; ++d) {
                for (int h = 0; h < 24; ++h) {
                    int k = d * 24 + h;
                    Rmx[k] = Rmxd[d]; Rmn[k] = Rmnd[d];
                    Rswmn[k] = Rswmnd[d]; Rlwmn[k] = Rlwmnd[d];
                    Rswmx[k] = Rswmxd[d]; Rlwmx[k] = Rlwmxd[d];
                    if (std::abs(Rnet[k]) > Gmxd[d]) Gmxd[d] = std::abs(Rnet[k]);
                }
                for (int h = 0; h < 24; ++h) Gmx[d * 24 + h] = Gmxd[d];
            }
            obspoint obstimeo{}; climpoint climo{}; vegpoint vegpo{};
            otherpoint othero{}; snowpoint snowo{};
            vegpo.pai = pai_m(i,j); vegpo.hgt = hgt_m(i,j); vegpo.clump = clump_m(i,j); vegpo.ltra = ltra_m(i,j);
            othero.slope = slope_m(i,j); othero.aspect = aspect_m(i,j);
            othero.lat = lats_m(i,j); othero.lon = lons_m(i,j);
            othero.zref = zref; othero.psim = 0.0; othero.psih = 0.0;
            int snowagec = isnowac_m(i, j);
            int snowageg = isnowag_m(i, j);
            double sdencp = ((sdp_r[0] - sdp_r[1]) * (1 - std::exp(-sdp_r[2] * isnowdc_m(i,j) / 100.0 -
                sdp_r[3] * snowagec / 24.0)) + sdp_r[1]) * 1000.0;
            double sdengp = ((sdp_r[0] - sdp_r[1]) * (1 - std::exp(-sdp_r[2] * isnowdg_m(i,j) * 0.5 / 100.0 -
                sdp_r[3] * snowageg / 24.0)) + sdp_r[1]) * 1000.0;
            double sdepcp = isnowdc_m(i,j);
            double sdepgp = isnowdg_m(i,j);
            meltc_w(i,j) = 0.0;
            meltg_w(i,j) = 0.0;
            for (int k = 0; k < tsteps; ++k) {
                int idx = i + rows * j + cols * rows * k;
                int snowtest = 0;
                if (sdepcp > 0.0) snowtest = 1;
                if (tc_v[idx] < 2.0 && prec_v[idx] > 0.0) snowtest = 1;
                if (snowtest > 0) {
                    double paip = pai_m(i,j);
                    if (hgt_m(i,j) > sdepgp) paip = paip * (hgt_m(i,j) - sdepgp) / hgt_m(i,j);
                    double dtR = Rmx[k] - Rmn[k];
                    double trS = skyview_m(i,j) * std::exp(-paip);
                    double Rem = 0.97 * sb * radem(tc_v[idx]);
                    double dmxS = trS * Rswmx[k] + trS * Rlwmx[k] + (1 - trS) * Rem - Rem;
                    double dmnS = trS * Rswmn[k] + trS * Rlwmn[k] + (1 - trS) * Rem - Rem;
                    double Gmu = (dmxS - dmnS) / dtR;
                    double G = Gp_v[idx] * Gmu;
                    if (G > Gmx[k]) G = Gmx[k];
                    if (G < -Gmx[k]) G = -Gmx[k];
                    solmodel solp = solpositionCpp(lats_m(i,j), lons_m(i,j), year_r[k], month_r[k], day_r[k], hour_r[k]);
                    int sindx = static_cast<int>(std::round(solp.azid / 15.0)) % 24;
                    double ha = hor_v[sindx * rows * cols + j * rows + i];
                    double sa = pi / 2.0 - solp.zenr;
                    double smu = 1.0;
                    if (ha > std::tan(sa)) smu = 0.0;
                    double ws = wsa_v[windex_r[k] * rows * cols + j * rows + i];
                    double u2p = umu_v[idx] * ws * u2_v[idx];
                    double Rdifp = Rdif_v[idx] * skyview_m(i,j);
                    double Rdirp = (Rsw_v[idx] - Rdif_v[idx]) * smu;
                    double Rswp = Rdirp + Rdifp;
                    double Rlwp = Rlw_v[idx] * skyview_m(i,j);
                    double ea = satvapCpp(tc_v[idx]) * rh_v[idx] / 100.0;
                    double te = (Tcp_v[idx] + tc_v[idx]) / 2.0;
                    obstimeo.year = year_r[k]; obstimeo.month = month_r[k];
                    obstimeo.day = day_r[k];   obstimeo.hour = hour_r[k];
                    climo.tc = tc_v[idx]; climo.ea = ea; climo.pk = pk_v[idx];
                    climo.u2 = u2p; climo.prec = prec_v[idx];
                    climo.Rsw = Rswp; climo.Rdif = Rdifp; climo.Rlw = Rlwp;
                    climo.Tci = Tcp_v[idx]; climo.te = te;
                    othero.G = G;
                    snowo.alb = salb[k]; snowo.sdenc = sdencp; snowo.sdeng = sdengp;
                    snowo.sdepc = sdepcp; snowo.sdepg = sdepgp;
                    snowo.snowagec = snowagec; snowo.snowageg = snowageg;
                    snowmodpoint smod = snowoneB(obstimeo, climo, vegpo, snowo, othero, sdp_r, 1.0);
                    Tc_w[idx] = smod.Tc; Tg_w[idx] = smod.Tg;
                    sdepc_w[idx] = smod.sdepc; sdepg_w[idx] = smod.sdepg; sdenc_w[idx] = smod.sdenc;
                    sdencp = smod.sdenc; sdengp = smod.sdeng;
                    sdepcp = smod.sdepc; sdepgp = smod.sdepg;
                    snowagec = smod.snowagec; snowageg = smod.snowageg;
                    double melc = smod.mSc + smod.mMc + smod.mRc;
                    double melg = smod.mSg + smod.mMg + smod.mRg;
                    meltc_w(i,j) = meltc_w(i,j) + (melc * 1000.0) / smod.sdenc;
                    meltg_w(i,j) = meltg_w(i,j) + (melg * 1000.0) / smod.sdeng;
                } else {
                    Tc_w[idx] = 0.0; Tg_w[idx] = 0.0;
                    sdepc_w[idx] = 0.0; sdepg_w[idx] = 0.0;
                    sdenc_w[idx] = sdp_r[1] * 1000.0;
                }
            }
            agec_w(i,j) = static_cast<double>(snowagec);
            ageg_w(i,j) = static_cast<double>(snowageg);
        }
    }
};

// [[Rcpp::export]]
List gridmodelsnow2Par(DataFrame obstime, List climdata, List pointm, List vegp,
    List other, std::string snowenv, int ncores = 2)
{
    IntegerVector year = obstime["year"]; IntegerVector month = obstime["month"];
    IntegerVector day = obstime["day"];   NumericVector hour = obstime["hour"];
    NumericMatrix pai = vegp["pai"]; NumericMatrix hgt = vegp["hgt"];
    NumericMatrix ltra = vegp["leaft"]; NumericMatrix clump = vegp["clump"];
    int rows = pai.nrow(); int cols = pai.ncol();
    int tsteps = year.size();
    IntegerVector dim = {rows, cols, tsteps};
    NumericVector Gp = pointm["Gp"]; NumericVector Tcp = pointm["Tc"];
    NumericVector RswabsG = pointm["RswabsG"]; NumericVector RlwabsG = pointm["RlwabsG"];
    NumericVector umu = pointm["umu"];
    NumericMatrix slope = other["slope"]; NumericMatrix aspect = other["aspect"];
    NumericMatrix skyview = other["skyview"];
    NumericVector wsa = other["wsa"]; NumericVector hor = other["hor"];
    NumericMatrix lats = other["lats"]; NumericMatrix lons = other["lons"];
    double zref = other["zref"];
    NumericMatrix isnowdc = other["isnowdc"]; NumericMatrix isnowdg = other["isnowdg"];
    IntegerMatrix isnowac = other["isnowac"]; IntegerMatrix isnowag = other["isnowag"];
    NumericVector tc = climdata["temp"]; NumericVector rh = climdata["relhum"];
    NumericVector pk = climdata["pres"]; NumericVector Rsw = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"]; NumericVector Rlw = climdata["lwdown"];
    NumericVector u2 = climdata["windspeed"]; NumericVector wdir = climdata["winddir"];
    NumericVector prec = climdata["precip"];
    std::vector<int> windex_v(tsteps);
    for (int k = 0; k < tsteps; ++k)
        windex_v[k] = static_cast<int>(std::round(wdir[k] / 45)) % 8;
    std::vector<double> sdp = snowdenp(snowenv);
    std::vector<int> yv(tsteps), mv(tsteps), dv(tsteps); std::vector<double> hv(tsteps);
    for (int k = 0; k < tsteps; ++k) { yv[k]=year[k]; mv[k]=month[k]; dv[k]=day[k]; hv[k]=hour[k]; }
    int n3 = rows * cols * tsteps;
    NumericVector Tc(n3, NA_REAL), Tg(n3, NA_REAL), sdepc(n3, NA_REAL), sdepg(n3, NA_REAL), sdenc(n3, NA_REAL);
    NumericMatrix agec = bioclimfill(rows, cols), ageg = bioclimfill(rows, cols);
    NumericMatrix meltc = bioclimfill(rows, cols), meltg = bioclimfill(rows, cols);
    Tc.attr("dim") = dim; Tg.attr("dim") = dim;
    sdepc.attr("dim") = dim; sdepg.attr("dim") = dim; sdenc.attr("dim") = dim;
    Snow2Worker worker(
        pai, hgt, ltra, clump, slope, aspect, skyview, lats, lons,
        isnowdc, isnowdg, isnowac, isnowag,
        wsa, hor, tc, rh, pk, Rsw, Rdif, Rlw, u2, prec,
        Gp, Tcp, RswabsG, RlwabsG, umu,
        windex_v, yv, mv, dv, hv, sdp,
        zref, rows, cols, tsteps,
        Tc, Tg, sdepc, sdepg, sdenc, agec, ageg, meltc, meltg);
    parallelFor(0, rows * cols, worker, 1, ncores);
    List out;
    out["Tc"] = Tc; out["Tg"] = Tg; out["sdepc"] = sdepc; out["sdepg"] = sdepg;
    out["sden"] = sdenc; out["agec"] = agec; out["ageg"] = ageg;
    out["meltc"] = meltc; out["meltg"] = meltg;
    return out;
}

// ============================================================
// GridMicroSnow1Worker: snow microclimate, data.frame climate
// ============================================================
struct GridMicroSnow1Worker : public Worker {
    const RMatrix<double> pai_m, paia_m, hgt_m, ltra_m, clump_m, leafd_m, leafden_m;
    const RMatrix<double> slope_m, aspect_m, skyview_m, Smax_m;
    const RVector<double> wsa_v, hor_v;
    const RVector<double> tc_v, rh_v, pk_v, Rsw_v, Rdif_v, Rlw_v, u2_v, umu_v;
    const RVector<double> snowtempc_v, snowtempg_v, swe_v, sdepg_v, sden_v, Tzd_v;
    const RMatrix<double> meanD_m;
    const RVector<double> salb_v;
    std::vector<double> zend_r, zenr_r, azid_r, azir_r;
    std::vector<int>    sindex_r, windex_r;
    double reqhgt, zref, mxtc, mat;
    bool Dynreqhgt;
    int rows, cols, tsteps, hiy;
    std::vector<bool> out_r;
    RVector<double> Tz_w, tleaf_w, relhum_w, soilm_w, uz_w;
    RVector<double> Rdirdown_w, Rdifdown_w, Rlwdown_w, Rswup_w, Rlwup_w;

    GridMicroSnow1Worker(
        NumericMatrix pai, NumericMatrix paia, NumericMatrix hgt,
        NumericMatrix ltra, NumericMatrix clump, NumericMatrix leafd, NumericMatrix leafden,
        NumericMatrix slope, NumericMatrix aspect, NumericMatrix skyview, NumericMatrix Smax,
        NumericVector wsa, NumericVector hor,
        NumericVector tc, NumericVector rh, NumericVector pk,
        NumericVector Rsw, NumericVector Rdif, NumericVector Rlw,
        NumericVector u2, NumericVector umu,
        NumericVector snowtempc, NumericVector snowtempg, NumericVector swe,
        NumericVector sdepg, NumericVector sden, NumericVector Tzd,
        NumericMatrix meanD, NumericVector salb,
        std::vector<double> zend, std::vector<double> zenr,
        std::vector<double> azid, std::vector<double> azir,
        std::vector<int> sindex, std::vector<int> windex,
        double reqhgt_, double zref_, double mxtc_, double mat_,
        bool Dynreqhgt_, int rows_, int cols_, int tsteps_, int hiy_,
        std::vector<bool> out,
        NumericVector Tz, NumericVector tleaf, NumericVector relhum, NumericVector soilm,
        NumericVector uz, NumericVector Rdirdown, NumericVector Rdifdown,
        NumericVector Rlwdown, NumericVector Rswup, NumericVector Rlwup
    ) :
        pai_m(pai), paia_m(paia), hgt_m(hgt), ltra_m(ltra), clump_m(clump),
        leafd_m(leafd), leafden_m(leafden),
        slope_m(slope), aspect_m(aspect), skyview_m(skyview), Smax_m(Smax),
        wsa_v(wsa), hor_v(hor),
        tc_v(tc), rh_v(rh), pk_v(pk), Rsw_v(Rsw), Rdif_v(Rdif), Rlw_v(Rlw),
        u2_v(u2), umu_v(umu),
        snowtempc_v(snowtempc), snowtempg_v(snowtempg), swe_v(swe),
        sdepg_v(sdepg), sden_v(sden), Tzd_v(Tzd),
        meanD_m(meanD), salb_v(salb),
        zend_r(std::move(zend)), zenr_r(std::move(zenr)),
        azid_r(std::move(azid)), azir_r(std::move(azir)),
        sindex_r(std::move(sindex)), windex_r(std::move(windex)),
        reqhgt(reqhgt_), zref(zref_), mxtc(mxtc_), mat(mat_),
        Dynreqhgt(Dynreqhgt_), rows(rows_), cols(cols_), tsteps(tsteps_), hiy(hiy_),
        out_r(std::move(out)),
        Tz_w(Tz), tleaf_w(tleaf), relhum_w(relhum), soilm_w(soilm), uz_w(uz),
        Rdirdown_w(Rdirdown), Rdifdown_w(Rdifdown), Rlwdown_w(Rlwdown),
        Rswup_w(Rswup), Rlwup_w(Rlwup)
    {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t cell = begin; cell < end; ++cell) {
            int i = static_cast<int>(cell) % rows;
            int j = static_cast<int>(cell) / rows;
            if (std::isnan(hgt_m(i, j))) continue;
            for (int k = 0; k < tsteps; ++k) {
                int idx = i + rows * j + cols * rows * k;
                if (swe_v[idx] > 0.0) {
                    double reqhgts = Dynreqhgt ? reqhgt : (reqhgt - sdepg_v[idx]);
                    if (reqhgts >= 0.0) {
                        int shadowmask = 1;
                        double ha = hor_v[sindex_r[k] * rows * cols + j * rows + i];
                        double sa = (pi / 2.0) - zenr_r[k];
                        double si = solarindexCpp(slope_m(i,j), aspect_m(i,j), zend_r[k], azid_r[k], true);
                        if (std::isnan(si)) si = std::cos(zend_r[k] * torad);
                        if (ha > std::tan(sa)) shadowmask = 0;
                        double ws = wsa_v[windex_r[k] * rows * cols + j * rows + i];
                        solmodel solp;
                        solp.zend = zend_r[k]; solp.zenr = zenr_r[k];
                        solp.azid = azid_r[k]; solp.azir = azir_r[k];
                        double sdepc = swe_v[idx] / sden_v[idx];
                        snowpoint2 snowp;
                        snowp.snowtempg = snowtempg_v[idx]; snowp.snowtempc = snowtempc_v[idx];
                        snowp.sdepc = sdepc; snowp.sdepg = sdepg_v[idx]; snowp.sdenc = sden_v[idx];
                        snowp.albc = salb_v[k]; snowp.albg = salb_v[k];
                        snowmicro apv = snowabovepoint(reqhgts, zref,
                            tc_v[k], rh_v[k], pk_v[k], u2_v[k],
                            Rsw_v[k], Rdif_v[k], Rlw_v[k],
                            hgt_m(i,j), pai_m(i,j), paia_m(i,j), leafd_m(i,j),
                            clump_m(i,j), ltra_m(i,j), leafden_m(i,j),
                            solp, si, skyview_m(i,j), shadowmask, ws, umu_v[k], mxtc, snowp);
                        if (out_r[0]) Tz_w[idx]        = apv.Tz;
                        if (out_r[1]) tleaf_w[idx]     = apv.tleaf;
                        if (out_r[2]) relhum_w[idx]    = apv.rh;
                        if (out_r[4]) uz_w[idx]        = apv.uz;
                        if (out_r[5]) Rdirdown_w[idx]  = apv.Rbdown;
                        if (out_r[6]) Rdifdown_w[idx]  = apv.Rddown;
                        if (out_r[7]) Rlwdown_w[idx]   = apv.Rlwdn;
                        if (out_r[8]) Rswup_w[idx]     = apv.Rdup;
                        if (out_r[9]) Rlwup_w[idx]     = apv.Rlwup;
                    } else {
                        double bpv = belowpointsnow(reqhgts, meanD_m(i,j),
                            snowtempg_v[idx], Tzd_v[idx], mat, hiy);
                        if (out_r[0]) Tz_w[idx]        = bpv;
                        if (out_r[1]) tleaf_w[idx]     = bpv;
                        if (out_r[2]) relhum_w[idx]    = 100.0;
                        if (out_r[4]) uz_w[idx]        = 0.0;
                        if (out_r[5]) Rdirdown_w[idx]  = 0.0;
                        if (out_r[6]) Rdifdown_w[idx]  = 0.0;
                        if (out_r[7]) Rlwdown_w[idx]   = 0.0;
                        if (out_r[8]) Rswup_w[idx]     = 0.0;
                        if (out_r[9]) Rlwup_w[idx]     = 0.0;
                    }
                    if (out_r[3]) soilm_w[idx] = Smax_m(i,j);
                }
            }
        }
    }
};

// [[Rcpp::export]]
List gridmicrosnow1Par(double reqhgt, bool Dynreqhgt, DataFrame obstime,
    DataFrame climdata, List snowm, List micro, List vegp, List other,
    double mat, std::vector<bool> out, int ncores = 2)
{
    IntegerVector year = obstime["year"]; IntegerVector month = obstime["month"];
    IntegerVector day  = obstime["day"];  NumericVector hour  = obstime["hour"];
    NumericVector tc   = climdata["temp"];    NumericVector rh   = climdata["relhum"];
    NumericVector pk   = climdata["pres"];    NumericVector Rsw  = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"];  NumericVector Rlw  = climdata["lwdown"];
    NumericVector u2   = climdata["windspeed"]; NumericVector wdir = climdata["winddir"];
    NumericVector prec = climdata["precip"]; NumericVector umu  = climdata["umu"];
    NumericMatrix pai     = vegp["pai"];   NumericMatrix paia    = vegp["paia"];
    NumericMatrix hgt     = vegp["hgt"];   NumericMatrix ltra    = vegp["leaft"];
    NumericMatrix clump   = vegp["clump"]; NumericMatrix leafd   = vegp["leafd"];
    NumericMatrix leafden = vegp["leafden"];
    NumericMatrix slope   = other["slope"];   NumericMatrix aspect  = other["aspect"];
    NumericMatrix skyview = other["skyview"]; NumericVector wsa     = other["wsa"];
    NumericVector hor     = other["hor"];     double lat   = other["lat"];
    double lon = other["lon"];               double zref  = other["zref"];
    NumericMatrix Smax    = other["Smax"];
    NumericVector snowtempc = snowm["Tc"];  NumericVector snowtempg = snowm["Tg"];
    NumericVector swe       = snowm["totalSWE"];
    NumericVector sdepg     = snowm["groundsnowdepth"];
    NumericVector sden      = snowm["snowden"];
    NumericVector Tzd  = snowdayan(snowm["Tg"]);
    NumericMatrix meanD = meanDsnow(snowm["snowden"]);
    int rows = pai.nrow(); int cols = pai.ncol(); int tsteps = tc.size();
    double mxtc = -273.15;
    std::vector<double> zend_v(tsteps), zenr_v(tsteps), azid_v(tsteps), azir_v(tsteps);
    std::vector<int>    sindex_v(tsteps), windex_v(tsteps);
    for (int k = 0; k < tsteps; ++k) {
        if (tc[k] > mxtc) mxtc = tc[k];
        solmodel sp = solpositionCpp(lat, lon, year[k], month[k], day[k], hour[k]);
        zend_v[k] = sp.zend; zenr_v[k] = sp.zenr;
        azid_v[k] = sp.azid; azir_v[k] = sp.azir;
        sindex_v[k] = static_cast<int>(std::round(sp.azid / 15)) % 24;
        windex_v[k] = static_cast<int>(std::round(wdir[k] / 45)) % 8;
    }
    NumericVector salb = snowalbCpp(prec);
    int hiy = (year[0] % 4 == 0 && (year[0] % 100 != 0 || year[0] % 400 == 0)) ? 366 * 24 : 365 * 24;
    int n3 = rows * cols * tsteps;
    IntegerVector dim3 = {rows, cols, tsteps};
    // Only allocate + seed the requested output layers from `micro`. Worker
    // writes are guarded by out_r[k]; unrequested outputs stay empty and are
    // never indexed. Avoids allocating 10 full arrays plus a copy of the input
    // when (typically) only out[c(1,4)] is requested in the snow path.
    NumericVector Tz_o, tleaf_o, relhum_o, soilm_o, uz_o, Rdir_o, Rdif_o, Rlw_o, Rswup_o, Rlwup_o;
    if (out[0]) { Tz_o     = NumericVector(n3, NA_REAL); NumericVector t = micro["Tz"];        std::copy(t.begin(), t.end(), Tz_o.begin()); }
    if (out[1]) { tleaf_o  = NumericVector(n3, NA_REAL); NumericVector t = micro["tleaf"];     std::copy(t.begin(), t.end(), tleaf_o.begin()); }
    if (out[2]) { relhum_o = NumericVector(n3, NA_REAL); NumericVector t = micro["relhum"];    std::copy(t.begin(), t.end(), relhum_o.begin()); }
    if (out[3]) { soilm_o  = NumericVector(n3, NA_REAL); NumericVector t = micro["soilm"];     std::copy(t.begin(), t.end(), soilm_o.begin()); }
    if (out[4]) { uz_o     = NumericVector(n3, NA_REAL); NumericVector t = micro["windspeed"]; std::copy(t.begin(), t.end(), uz_o.begin()); }
    if (out[5]) { Rdir_o   = NumericVector(n3, NA_REAL); NumericVector t = micro["Rdirdown"];  std::copy(t.begin(), t.end(), Rdir_o.begin()); }
    if (out[6]) { Rdif_o   = NumericVector(n3, NA_REAL); NumericVector t = micro["Rdifdown"];  std::copy(t.begin(), t.end(), Rdif_o.begin()); }
    if (out[7]) { Rlw_o    = NumericVector(n3, NA_REAL); NumericVector t = micro["Rlwdown"];   std::copy(t.begin(), t.end(), Rlw_o.begin()); }
    if (out[8]) { Rswup_o  = NumericVector(n3, NA_REAL); NumericVector t = micro["Rswup"];     std::copy(t.begin(), t.end(), Rswup_o.begin()); }
    if (out[9]) { Rlwup_o  = NumericVector(n3, NA_REAL); NumericVector t = micro["Rlwup"];     std::copy(t.begin(), t.end(), Rlwup_o.begin()); }
    GridMicroSnow1Worker worker(
        pai, paia, hgt, ltra, clump, leafd, leafden,
        slope, aspect, skyview, Smax, wsa, hor,
        tc, rh, pk, Rsw, Rdif, Rlw, u2, umu,
        snowtempc, snowtempg, swe, sdepg, sden, Tzd, meanD, salb,
        zend_v, zenr_v, azid_v, azir_v, sindex_v, windex_v,
        reqhgt, zref, mxtc, mat, Dynreqhgt, rows, cols, tsteps, hiy, out,
        Tz_o, tleaf_o, relhum_o, soilm_o, uz_o,
        Rdir_o, Rdif_o, Rlw_o, Rswup_o, Rlwup_o);
    parallelFor(0, rows * cols, worker, 1, ncores);
    List outp;
    if (out[0]) { Tz_o.attr("dim")     = dim3; outp["Tz"]        = Tz_o; }
    if (out[1]) { tleaf_o.attr("dim")  = dim3; outp["tleaf"]     = tleaf_o; }
    if (out[2]) { relhum_o.attr("dim") = dim3; outp["relhum"]    = relhum_o; }
    if (out[3]) { soilm_o.attr("dim")  = dim3; outp["soilm"]     = soilm_o; }
    if (out[4]) { uz_o.attr("dim")     = dim3; outp["windspeed"] = uz_o; }
    if (out[5]) { Rdir_o.attr("dim")   = dim3; outp["Rdirdown"]  = Rdir_o; }
    if (out[6]) { Rdif_o.attr("dim")   = dim3; outp["Rdifdown"]  = Rdif_o; }
    if (out[7]) { Rlw_o.attr("dim")    = dim3; outp["Rlwdown"]   = Rlw_o; }
    if (out[8]) { Rswup_o.attr("dim")  = dim3; outp["Rswup"]     = Rswup_o; }
    if (out[9]) { Rlwup_o.attr("dim")  = dim3; outp["Rlwup"]     = Rlwup_o; }
    return outp;
}

// ============================================================
// GridMicroSnow2Worker: snow microclimate, array climate
// ============================================================
struct GridMicroSnow2Worker : public Worker {
    const RMatrix<double> pai_m, paia_m, hgt_m, ltra_m, clump_m, leafd_m, leafden_m;
    const RMatrix<double> slope_m, aspect_m, skyview_m, lats_m, lons_m, Smax_m;
    const RVector<double> wsa_v, hor_v;
    const RVector<double> tc_v, rh_v, pk_v, Rsw_v, Rdif_v, Rlw_v, u2_v, umu_v, prec_v;
    const RVector<double> snowtempc_v, snowtempg_v, swe_v, sdepg_v, sden_v, Tzd_v;
    const RMatrix<double> meanD_m;
    std::vector<int>    windex_r, year_r, month_r, day_r;
    std::vector<double> hour_r;
    double reqhgt, zref, mat;
    bool Dynreqhgt;
    int rows, cols, tsteps, hiy;
    std::vector<bool> out_r;
    RVector<double> Tz_w, tleaf_w, relhum_w, soilm_w, uz_w;
    RVector<double> Rdirdown_w, Rdifdown_w, Rlwdown_w, Rswup_w, Rlwup_w;

    GridMicroSnow2Worker(
        NumericMatrix pai, NumericMatrix paia, NumericMatrix hgt,
        NumericMatrix ltra, NumericMatrix clump, NumericMatrix leafd, NumericMatrix leafden,
        NumericMatrix slope, NumericMatrix aspect, NumericMatrix skyview,
        NumericMatrix lats, NumericMatrix lons, NumericMatrix Smax,
        NumericVector wsa, NumericVector hor,
        NumericVector tc, NumericVector rh, NumericVector pk,
        NumericVector Rsw, NumericVector Rdif, NumericVector Rlw,
        NumericVector u2, NumericVector umu, NumericVector prec,
        NumericVector snowtempc, NumericVector snowtempg, NumericVector swe,
        NumericVector sdepg, NumericVector sden, NumericVector Tzd,
        NumericMatrix meanD,
        std::vector<int> windex, std::vector<int> year, std::vector<int> month,
        std::vector<int> day, std::vector<double> hour,
        double reqhgt_, double zref_, double mat_,
        bool Dynreqhgt_, int rows_, int cols_, int tsteps_, int hiy_,
        std::vector<bool> out,
        NumericVector Tz, NumericVector tleaf, NumericVector relhum, NumericVector soilm,
        NumericVector uz, NumericVector Rdirdown, NumericVector Rdifdown,
        NumericVector Rlwdown, NumericVector Rswup, NumericVector Rlwup
    ) :
        pai_m(pai), paia_m(paia), hgt_m(hgt), ltra_m(ltra), clump_m(clump),
        leafd_m(leafd), leafden_m(leafden),
        slope_m(slope), aspect_m(aspect), skyview_m(skyview),
        lats_m(lats), lons_m(lons), Smax_m(Smax),
        wsa_v(wsa), hor_v(hor),
        tc_v(tc), rh_v(rh), pk_v(pk), Rsw_v(Rsw), Rdif_v(Rdif), Rlw_v(Rlw),
        u2_v(u2), umu_v(umu), prec_v(prec),
        snowtempc_v(snowtempc), snowtempg_v(snowtempg), swe_v(swe),
        sdepg_v(sdepg), sden_v(sden), Tzd_v(Tzd),
        meanD_m(meanD),
        windex_r(std::move(windex)), year_r(std::move(year)), month_r(std::move(month)),
        day_r(std::move(day)), hour_r(std::move(hour)),
        reqhgt(reqhgt_), zref(zref_), mat(mat_),
        Dynreqhgt(Dynreqhgt_), rows(rows_), cols(cols_), tsteps(tsteps_), hiy(hiy_),
        out_r(std::move(out)),
        Tz_w(Tz), tleaf_w(tleaf), relhum_w(relhum), soilm_w(soilm), uz_w(uz),
        Rdirdown_w(Rdirdown), Rdifdown_w(Rdifdown), Rlwdown_w(Rlwdown),
        Rswup_w(Rswup), Rlwup_w(Rlwup)
    {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t cell = begin; cell < end; ++cell) {
            int i = static_cast<int>(cell) % rows;
            int j = static_cast<int>(cell) / rows;
            if (std::isnan(hgt_m(i, j))) continue;
            double mxtc = -273.15;
            std::vector<int> hs(tsteps, 0);
            for (int k = 0; k < tsteps; ++k) {
                int idx = i + rows * j + cols * rows * k;
                if (tc_v[idx] > mxtc) mxtc = tc_v[idx];
                hs[k] = (k > 0) ? ((prec_v[idx] > 0.0) ? 0 : hs[k-1] + 1) : 0;
            }
            std::vector<double> salb(tsteps);
            for (int k = 0; k < tsteps; ++k) {
                // integer division hs[k]/24 matches snowalbCpp: log(0)=-inf → clamped to 0.95 for first 24h
                salb[k] = (-9.8740 * std::log(static_cast<double>(hs[k] / 24)) + 78.3434) / 100.0;
                if (salb[k] > 0.95) salb[k] = 0.95;
                if (salb[k] < 0.1)  salb[k] = 0.1;
            }
            for (int k = 0; k < tsteps; ++k) {
                int idx = i + rows * j + cols * rows * k;
                if (swe_v[idx] > 0.0) {
                    double reqhgts = Dynreqhgt ? reqhgt : (reqhgt - sdepg_v[idx]);
                    if (reqhgts >= 0.0) {
                        solmodel solp = solpositionCpp(lats_m(i,j), lons_m(i,j),
                            year_r[k], month_r[k], day_r[k], hour_r[k]);
                        int sindx = static_cast<int>(std::round(solp.azid / 15)) % 24;
                        int shadowmask = 1;
                        double ha = hor_v[sindx * rows * cols + j * rows + i];
                        double sa = (pi / 2.0) - solp.zenr;
                        double si = solarindexCpp(slope_m(i,j), aspect_m(i,j), solp.zend, solp.azid, true);
                        if (std::isnan(si)) si = std::cos(solp.zenr);
                        if (ha > std::tan(sa)) shadowmask = 0;
                        double ws = wsa_v[windex_r[k] * rows * cols + j * rows + i];
                        double sdepc = swe_v[idx] / sden_v[idx];
                        snowpoint2 snowp;
                        snowp.snowtempg = snowtempg_v[idx]; snowp.snowtempc = snowtempc_v[idx];
                        snowp.sdepc = sdepc; snowp.sdepg = sdepg_v[idx]; snowp.sdenc = sden_v[idx];
                        snowp.albc = salb[k]; snowp.albg = salb[k];
                        snowmicro apv = snowabovepoint(reqhgts, zref,
                            tc_v[idx], rh_v[idx], pk_v[idx], u2_v[idx],
                            Rsw_v[idx], Rdif_v[idx], Rlw_v[idx],
                            hgt_m(i,j), pai_m(i,j), paia_m(i,j), leafd_m(i,j),
                            clump_m(i,j), ltra_m(i,j), leafden_m(i,j),
                            solp, si, skyview_m(i,j), shadowmask, ws, umu_v[idx], mxtc, snowp);
                        if (out_r[0]) Tz_w[idx]        = apv.Tz;
                        if (out_r[1]) tleaf_w[idx]     = apv.tleaf;
                        if (out_r[2]) relhum_w[idx]    = apv.rh;
                        if (out_r[4]) uz_w[idx]        = apv.uz;
                        if (out_r[5]) Rdirdown_w[idx]  = apv.Rbdown;
                        if (out_r[6]) Rdifdown_w[idx]  = apv.Rddown;
                        if (out_r[7]) Rlwdown_w[idx]   = apv.Rlwdn;
                        if (out_r[8]) Rswup_w[idx]     = apv.Rdup;
                        if (out_r[9]) Rlwup_w[idx]     = apv.Rlwup;
                    } else {
                        double bpv = belowpointsnow(reqhgts, meanD_m(i,j),
                            snowtempg_v[idx], Tzd_v[idx], mat, hiy);
                        if (out_r[0]) Tz_w[idx]        = bpv;
                        if (out_r[1]) tleaf_w[idx]     = bpv;
                        if (out_r[2]) relhum_w[idx]    = 100.0;
                        if (out_r[4]) uz_w[idx]        = 0.0;
                        if (out_r[5]) Rdirdown_w[idx]  = 0.0;
                        if (out_r[6]) Rdifdown_w[idx]  = 0.0;
                        if (out_r[7]) Rlwdown_w[idx]   = 0.0;
                        if (out_r[8]) Rswup_w[idx]     = 0.0;
                        if (out_r[9]) Rlwup_w[idx]     = 0.0;
                    }
                    if (out_r[3]) soilm_w[idx] = Smax_m(i,j);
                }
            }
        }
    }
};

// [[Rcpp::export]]
List gridmicrosnow2Par(double reqhgt, bool Dynreqhgt, DataFrame obstime,
    List climdata, List snowm, List micro, List vegp, List other,
    double mat, std::vector<bool> out, int ncores = 2)
{
    IntegerVector year = obstime["year"]; IntegerVector month = obstime["month"];
    IntegerVector day  = obstime["day"];  NumericVector hour  = obstime["hour"];
    NumericVector tc   = climdata["temp"];    NumericVector rh   = climdata["relhum"];
    NumericVector pk   = climdata["pres"];    NumericVector Rsw  = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"];  NumericVector Rlw  = climdata["lwdown"];
    NumericVector u2   = climdata["windspeed"]; NumericVector wdir = climdata["winddir"];
    NumericVector prec = climdata["prec"];   NumericVector umu  = climdata["umu"];
    NumericMatrix pai     = vegp["pai"];   NumericMatrix paia    = vegp["paia"];
    NumericMatrix hgt     = vegp["hgt"];   NumericMatrix ltra    = vegp["leaft"];
    NumericMatrix clump   = vegp["clump"]; NumericMatrix leafd   = vegp["leafd"];
    NumericMatrix leafden = vegp["leafden"];
    NumericMatrix slope   = other["slope"];   NumericMatrix aspect  = other["aspect"];
    NumericMatrix skyview = other["skyview"]; NumericVector wsa     = other["wsa"];
    NumericVector hor     = other["hor"];     NumericMatrix lats    = other["lat"];
    NumericMatrix lons    = other["lon"];     double zref           = other["zref"];
    NumericMatrix Smax    = other["Smax"];
    NumericVector snowtempc = snowm["Tc"];  NumericVector snowtempg = snowm["Tg"];
    NumericVector swe       = snowm["totalSWE"];
    NumericVector sdepg     = snowm["groundsnowdepth"];
    NumericVector sden      = snowm["snowden"];
    NumericVector Tzd  = snowdayan(snowm["Tg"]);
    NumericMatrix meanD = meanDsnow(snowm["snowden"]);
    int rows = pai.nrow(); int cols = pai.ncol(); int tsteps = hour.size();
    int hiy = (year[0] % 4 == 0 && (year[0] % 100 != 0 || year[0] % 400 == 0)) ? 366 * 24 : 365 * 24;
    std::vector<int> windex_v(tsteps);
    for (int k = 0; k < tsteps; ++k)
        windex_v[k] = static_cast<int>(std::round(wdir[k] / 45)) % 8;
    std::vector<int> yv(tsteps), mv(tsteps), dv(tsteps); std::vector<double> hv(tsteps);
    for (int k = 0; k < tsteps; ++k) { yv[k]=year[k]; mv[k]=month[k]; dv[k]=day[k]; hv[k]=hour[k]; }
    int n3 = rows * cols * tsteps;
    IntegerVector dim3 = {rows, cols, tsteps};
    // Only allocate + seed the requested output layers from `micro`. Worker
    // writes are guarded by out_r[k]; unrequested outputs stay empty and are
    // never indexed. Avoids allocating 10 full arrays plus a copy of the input
    // when (typically) only out[c(1,4)] is requested in the snow path.
    NumericVector Tz_o, tleaf_o, relhum_o, soilm_o, uz_o, Rdir_o, Rdif_o, Rlw_o, Rswup_o, Rlwup_o;
    if (out[0]) { Tz_o     = NumericVector(n3, NA_REAL); NumericVector t = micro["Tz"];        std::copy(t.begin(), t.end(), Tz_o.begin()); }
    if (out[1]) { tleaf_o  = NumericVector(n3, NA_REAL); NumericVector t = micro["tleaf"];     std::copy(t.begin(), t.end(), tleaf_o.begin()); }
    if (out[2]) { relhum_o = NumericVector(n3, NA_REAL); NumericVector t = micro["relhum"];    std::copy(t.begin(), t.end(), relhum_o.begin()); }
    if (out[3]) { soilm_o  = NumericVector(n3, NA_REAL); NumericVector t = micro["soilm"];     std::copy(t.begin(), t.end(), soilm_o.begin()); }
    if (out[4]) { uz_o     = NumericVector(n3, NA_REAL); NumericVector t = micro["windspeed"]; std::copy(t.begin(), t.end(), uz_o.begin()); }
    if (out[5]) { Rdir_o   = NumericVector(n3, NA_REAL); NumericVector t = micro["Rdirdown"];  std::copy(t.begin(), t.end(), Rdir_o.begin()); }
    if (out[6]) { Rdif_o   = NumericVector(n3, NA_REAL); NumericVector t = micro["Rdifdown"];  std::copy(t.begin(), t.end(), Rdif_o.begin()); }
    if (out[7]) { Rlw_o    = NumericVector(n3, NA_REAL); NumericVector t = micro["Rlwdown"];   std::copy(t.begin(), t.end(), Rlw_o.begin()); }
    if (out[8]) { Rswup_o  = NumericVector(n3, NA_REAL); NumericVector t = micro["Rswup"];     std::copy(t.begin(), t.end(), Rswup_o.begin()); }
    if (out[9]) { Rlwup_o  = NumericVector(n3, NA_REAL); NumericVector t = micro["Rlwup"];     std::copy(t.begin(), t.end(), Rlwup_o.begin()); }
    GridMicroSnow2Worker worker(
        pai, paia, hgt, ltra, clump, leafd, leafden,
        slope, aspect, skyview, lats, lons, Smax, wsa, hor,
        tc, rh, pk, Rsw, Rdif, Rlw, u2, umu, prec,
        snowtempc, snowtempg, swe, sdepg, sden, Tzd, meanD,
        windex_v, yv, mv, dv, hv,
        reqhgt, zref, mat, Dynreqhgt, rows, cols, tsteps, hiy, out,
        Tz_o, tleaf_o, relhum_o, soilm_o, uz_o,
        Rdir_o, Rdif_o, Rlw_o, Rswup_o, Rlwup_o);
    parallelFor(0, rows * cols, worker, 1, ncores);
    List outp;
    if (out[0]) { Tz_o.attr("dim")     = dim3; outp["Tz"]        = Tz_o; }
    if (out[1]) { tleaf_o.attr("dim")  = dim3; outp["tleaf"]     = tleaf_o; }
    if (out[2]) { relhum_o.attr("dim") = dim3; outp["relhum"]    = relhum_o; }
    if (out[3]) { soilm_o.attr("dim")  = dim3; outp["soilm"]     = soilm_o; }
    if (out[4]) { uz_o.attr("dim")     = dim3; outp["windspeed"] = uz_o; }
    if (out[5]) { Rdir_o.attr("dim")   = dim3; outp["Rdirdown"]  = Rdir_o; }
    if (out[6]) { Rdif_o.attr("dim")   = dim3; outp["Rdifdown"]  = Rdif_o; }
    if (out[7]) { Rlw_o.attr("dim")    = dim3; outp["Rlwdown"]   = Rlw_o; }
    if (out[8]) { Rswup_o.attr("dim")  = dim3; outp["Rswup"]     = Rswup_o; }
    if (out[9]) { Rlwup_o.attr("dim")  = dim3; outp["Rlwup"]     = Rlwup_o; }
    return outp;
}

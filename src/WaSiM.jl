module WaSiM

src_path = joinpath(@__DIR__, "src")
# Include submodules
# include("DataProcessing.jl")
# include("Utils.jl")
include("smallfuncs.jl")
include("func-win.jl")
include("rasterfuncs.jl")
# include("pyjl.jl")
# include("cairomakie.jl")
# include("rcall.jl")
# include("wajs.jl")
# using .wa
# using .rst
#using .pyjl
#using .cmk
# using .rr
# using .wajs

#
#
# include("ncremover.jl")
#include("wasim-fncs.jl")

# Using submodules
# using .DataProcessing
# using .Utils
# using .FuncWin
# using .WaSimFncs

# Export functions
# DATAFRAME Operations
export dfrib, dfsp, dfsplog, dfp, dfp!, dfpall,
    dfr,
    yrsum, yrmean, dfm, dfmo, dfl, dfl!, dfilter,
    monmean, monsum, lastcol
#dfroute,

# File Operations
export readall, read_log_file,
    read_soildata, read_soildata_2, read_soildata_4,
    read_until_flag, read_wq,
    pout,
    read_landuse_data2,
    readroute, sf

# Hydrological Functions
export hd, hydro, hydro_f, hydromon, hyeval,
    kge, kge1, kge2, kge_df, kge_df3, kge_fread, kge_read,
    kge_rec, kgedf, kgerec, kgeval,
    yearsum, hyd2

# Statistical Functions
export ls, du, addname, all_values_equal,
    barsum, baryr,
    baryrmean, baryrsum, build_soil_dictionary, byear, climateplot,
    climateplot_dfold, cmplot, cntcolread, cntcols, colorfunction, colsums,
    condasize, corrbar, correlogram, ddense, denselog, dfilter, dfl, dfl!, dfm, dfmo, dfp,
    dfsp, dfsplog, dftrend,
    dprbig, dsbar, dtrange, dubplot, eeread, extract_duration_from_xml,
    extract_layers, f1_score, fdf, fillmissing, filterdf, findalls, findctl, fparse, fread, freaddf, fsoil, generate_export_statement, getdf, getf, getnames,
    ggofbatch, ggofjl, ggofjl_nosvg, globdf, globf, gofbatch, gofbatch_nosvg,
    grec, grep_files, grep_in_files, grep_with_context, gwread, hd, heat,
    heatraw, homes, ht, hydro, hydro_f, hydromon, hyeval, irfan, isroute, jldf,
    jldfnm, jsrast, jsread, juliasize, kernelplot, kge, kge1, kge2, kge_df,
    kge_df3, kge_fread, kge_read, kge_rec, kgedf, kgerec, kgeval, klog, ldf,
    ldfpall, lg, listdfs, loadalldfs, loadso, lplot, lplotf, lpro,
    luheat, lutab, luvars, mall, malldf, mbx,
    moisture_plot_with_confidence, monc, monc_f, moncol, monp
#revcrds, reverse_coords #moved to rst

# Plotting Functions
export aplot, bardf, bardfm, bardfm!, bargroup, barp,
    dfrib, dfp!, dfpall, ftp, ftplin, ftsp, filterplot!, dpr, dpr!, filterplot, denseplot,
    dfbar, stplot!, moisture_plot_with_confidence, mbx, lplot, lplotf, lpro, luheat, luplot, luscatter

export lutab, luvars, mall, malldf,
    monc, monc_f, moncol, monmean, monp, monsum,
    mvwasim2, mywd, nctodf, nctodfo, nqp, nread, nse, nse2, nse_rec, nseval,
    nsevalraw, nsx, odfr, pall, penman_monteith, pers, pout, plot_duration_graph,
    plot_grouped_metrics, plotf, polygonize_raster, print_lines_between_patterns,
    process_file2, process_folders_and_subfolders, pw, pxm,
    qall_num, qba, qbb, qpl, qplot, qrtr, rds, read_landuse_data2,
    read_log_file, read_soildata, read_soildata_2, read_soildata_4,
    read_soildata_raw, read_until_flag, read_wq, readall, readdf, readf,
    readmhm, readroute, rec, rename_columns, rename_columns!, rename_duplicates,
    renamer, reorder_df, rmm, rmopt, rmqout, routeg, routg,
    rowmeans, rowsums, runwasim, seldf, selt, skipyr, so_read, stplot!, strans

end

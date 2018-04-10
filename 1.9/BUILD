package(default_visibility = ["//visibility:public"])

cc_binary(
    name = "plink",
    srcs = [
        "plink.c",
        "plink_assoc.c",
        "plink_calc.c",
        "plink_cluster.c",
        "plink_cnv.c",
        "plink_common.c",
        "plink_data.c",
        "plink_dosage.c",
        "plink_family.c",
        "plink_filter.c",
        "plink_glm.c",
        "plink_help.c",
        "plink_homozyg.c",
        "plink_lasso.c",
        "plink_ld.c",
        "plink_matrix.c",
        "plink_misc.c",
        "plink_perm.c",
        "plink_rserve.c",
        "plink_set.c",
        "plink_stats.c",
        "SFMT.c",
        "dcdflib.c",
        "pigz.c",
        "yarn.c",
        "hfile.c",
        "bgzf.c",
        "Rconnection.cc"
    ] + glob(["**/*.h"]),
    copts = [
        "-x c++",
        "-DDYNAMIC_ZLIB"
    ],
    linkopts = [
        "-framework Accelerate",
        "-lpthread",
        "-ldl",
        "-lz",
    ]
)
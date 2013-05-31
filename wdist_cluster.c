#include "wdist_cluster.h"

void cluster_init(Cluster_info* cluster_ptr) {
  cluster_ptr->fname = NULL;
  cluster_ptr->match_fname = NULL;
  cluster_ptr->match_type_fname = NULL;
  cluster_ptr->qmatch_fname = NULL;
  cluster_ptr->qt_fname = NULL;
  cluster_ptr->modifier = 0;
  cluster_ptr->ppc = 0.0;
  cluster_ptr->max_size = 0xffffffffU;
  cluster_ptr->max_cases = 0xffffffffU;
  cluster_ptr->max_controls = 0xffffffffU;
  cluster_ptr->min_ct = 1;
  cluster_ptr->mds_dim_ct = 0;
  cluster_ptr->neighbor_n1 = 0; // yeah, American spelling, deal with it
  cluster_ptr->neighbor_n2 = 0;
  cluster_ptr->max_missing_discordance = 1.0;
}

#include "workspace.h"

blsp_workspace *blsp_workspace_sampler_alloc(const size_t nobs,
                                             const size_t nyrs,
                                             const size_t nburn,
                                             const size_t niter,
                                             const size_t npar) {

  // Parameter checking
  if (niter <= nburn) {
    GSL_ERROR(
        "number of burn iterations must be less than number of iterations", 0);
  }
  if (nobs == 0) {
    GSL_ERROR("number of observations must be > 0", 0);
  }

  if (nyrs == 0) {
    GSL_ERROR("number of years must be > 0", 0);
  }

  // NOTE: Surely there is a way to hack together a memory arena based of the
  //       gsl_block structure. That would drastically reduce the number of
  //       memory allocs needed and would allow the workspace to be freed
  //       just by freeing the underlying memory block (at least I think)
  blsp_workspace *w;

  w = calloc(1, sizeof(blsp_workspace));

  if (w == 0) {
    GSL_ERROR_VAL("failed to allocate space for blsp struct", GSL_ENOMEM, 0);
  }

  w->nobs = nobs; /* number of observations */
  w->nyrs = nyrs; /* number of year */
  w->npar = npar; /* number of parameters */
  w->niter = niter;
  w->nburn = nburn;
  w->sigma = 0.08040725; /* global sttdev for all parameters in all years */

  w->theta_hat = gsl_matrix_alloc(nyrs, npar);

  if (w->theta_hat == 0) {
    blsp_free_workspace(w);
    GSL_ERROR_VAL("failed to allocate space for theta_hat matrix", GSL_ENOMEM,
                  0);
  }

  w->mh_mat = gsl_matrix_alloc(nyrs, npar);

  if (w->mh_mat == 0) {
    blsp_free_workspace(w);
    GSL_ERROR_VAL("failed to allocate space for mh_mat matrix", GSL_ENOMEM, 0);
  }

  w->accept_mat = gsl_matrix_alloc(nyrs, npar);

  if (w->accept_mat == 0) {
    blsp_free_workspace(w);
    GSL_ERROR_VAL("failed to allocate space for accept_mat matrix", GSL_ENOMEM,
                  0);
  }
  gsl_matrix_set_all(w->accept_mat, 0.1);

  w->attempt_mat = gsl_matrix_alloc(nyrs, npar);

  if (w->attempt_mat == 0) {
    blsp_free_workspace(w);
    GSL_ERROR_VAL("failed to allocate space for attempt_mat matrix", GSL_ENOMEM,
                  0);
  }
  gsl_matrix_set_all(w->attempt_mat, 0.1);

  w->parameter_track = gsl_matrix_alloc(nyrs * w->niter, 7 + 2);

  if (w->parameter_track == 0) {
    blsp_free_workspace(w);
    GSL_ERROR_VAL("failed to allocate space for parameter tracking matrix",
                  GSL_ENOMEM, 0);
  }

  w->vi_vec = gsl_vector_alloc(nobs);

  if (w->vi_vec == 0) {
    blsp_free_workspace(w);
    GSL_ERROR_VAL("failed to allocate space for vi_vec matrix", GSL_ENOMEM, 0);
  }

  w->doy_vec = gsl_vector_alloc(nobs);

  if (w->doy_vec == 0) {
    blsp_free_workspace(w);
    GSL_ERROR_VAL("failed to allocate space for doy_vec matrix", GSL_ENOMEM, 0);
  }

  w->can_p_vec = gsl_vector_alloc(npar);

  if (w->can_p_vec == 0) {
    blsp_free_workspace(w);
    GSL_ERROR_VAL("failed to allocate space for can_p_vec matrix", GSL_ENOMEM,
                  0);
  }

  w->theta_mu_vec = gsl_vector_calloc(npar);

  if (w->theta_mu_vec == 0) {
    blsp_free_workspace(w);
    GSL_ERROR_VAL("failed to allocate space for theta_mu_vec matrix",
                  GSL_ENOMEM, 0);
  }

  w->theta_sd_vec = gsl_vector_calloc(npar);

  if (w->theta_sd_vec == 0) {
    blsp_free_workspace(w);
    GSL_ERROR_VAL("failed to allocate space for theta_sd_vec matrix",
                  GSL_ENOMEM, 0);
  }

  w->loglike_vec = gsl_vector_alloc(nyrs);

  if (w->loglike_vec == 0) {
    blsp_free_workspace(w);
    GSL_ERROR_VAL("failed to allocate space for loglike_vec matrix", GSL_ENOMEM,
                  0);
  }

  return w;
}

void blsp_free_workspace(blsp_workspace *w) {

  RETURN_IF_NULL(w);

  if (w->theta_hat)
    gsl_matrix_free(w->theta_hat);

  if (w->mh_mat)
    gsl_matrix_free(w->mh_mat);

  if (w->accept_mat)
    gsl_matrix_free(w->accept_mat);

  if (w->attempt_mat)
    gsl_matrix_free(w->attempt_mat);

  if (w->parameter_track)
    gsl_matrix_free(w->parameter_track);

  if (w->vi_vec)
    gsl_vector_free(w->vi_vec);

  if (w->doy_vec)
    gsl_vector_free(w->doy_vec);

  if (w->can_p_vec)
    gsl_vector_free(w->can_p_vec);

  if (w->theta_mu_vec)
    gsl_vector_free(w->theta_mu_vec);

  if (w->theta_sd_vec)
    gsl_vector_free(w->theta_sd_vec);

  if (w->loglike_vec)
    gsl_vector_free(w->loglike_vec);

  free(w);
}

void fill_matrix_by_row(gsl_matrix *m, gsl_vector *x) {
  size_t nrow;
  nrow = m->size1;
  for (size_t i = 0; i < nrow; ++i)
    gsl_matrix_set_row(m, i, x);
}

void workspace_init_thetas(gsl_vector *mu, gsl_vector *sd, blsp_workspace *w) {
  gsl_vector_memcpy(w->theta_mu_vec, mu);
  gsl_vector_memcpy(w->theta_sd_vec, sd);
  fill_matrix_by_row(w->theta_hat, w->theta_mu_vec);
  fill_matrix_by_row(w->mh_mat, w->theta_sd_vec);
}

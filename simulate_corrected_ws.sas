/* Simulation study for the corrected Welch-Satterthwaite estimator.

   Run this program from the repository directory, or change the %include
   path below to point to corrected_welch_satterthwaite_iml.sas.

   The SAS random-number stream is not the same as the R or Python stream,
   so the Monte Carlo estimates are structurally comparable but not
   numerically identical to the output shown for the other languages.
*/

proc iml;
   %include "corrected_welch_satterthwaite_iml.sas";

   seed      = 20260715;
   K_values  = {2 4 6 8 10 20 30 40};
   nu_values = {1 2 3 4 5 10};
   nsim      = 10000;

   nres    = ncol(K_values) * ncol(nu_values);
   results = J(nres, 8, .);

   call randseed(seed);

   row = 1;
   do ik = 1 to ncol(K_values);
      K = K_values[ik];

      do inu = 1 to ncol(nu_values);
         nu     = nu_values[inu];
         target = K * nu;
         ones   = J(K, 1, 1);

         /* Draw nsim K-vectors of chi-square variance components. */
         x  = J(nsim, K, .);
         call randgen(x, 'CHISQUARE', nu);
         s2 = x #/ nu;

         /* Row-wise corrected Welch-Satterthwaite estimate. */
         num = s2 * ones;
         den = ((s2 # s2) #/ J(nsim, K, nu + 2)) * ones;
         df_est = (num # num) #/ den - 2;

         mean_df = sum(df_est) / nsim;
         bias    = mean_df - target;
         rel_bias = bias / target;

         dev     = df_est - J(nsim, 1, mean_df);
         mc_se   = sqrt(sum(dev # dev) / (nsim - 1)) / sqrt(nsim);

         dev_t   = df_est - J(nsim, 1, target);
         rmse    = sqrt(sum(dev_t # dev_t) / nsim);

         results[row, 1] = K;
         results[row, 2] = nu;
         results[row, 3] = target;
         results[row, 4] = mean_df;
         results[row, 5] = bias;
         results[row, 6] = rel_bias;
         results[row, 7] = mc_se;
         results[row, 8] = rmse;

         row = row + 1;
      end;
   end;

   print results
      [COLNAME={'K' 'nu' 'target_df' 'mean_df' 'bias' 'relative_bias'
                'mc_se' 'rmse'}];
quit;

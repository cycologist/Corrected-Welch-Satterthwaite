/* Corrected Welch-Satterthwaite SAS/IML modules.

   Include this file inside PROC IML:

     proc iml;
        %include "corrected_welch_satterthwaite_iml.sas";

        s2 = {2 3};
        nu = {10 20};
        w  = {1 2};

        df1 = cws2(2, 3, 10, 20, 1, 1);
        df2 = cws(s2, nu, w);

        print df1 df2;
     quit;

   The corrected estimator is:

     (sum_k w_k S_k^2)^2 / sum_k w_k^2 S_k^4 / (nu_k + 2) - 2
*/

start cws(s2, nu, w);
   s2 = colvec(s2);
   nu = colvec(nu);
   w  = colvec(w);
   n  = nrow(s2);

   if n <= 0 then do;
      print "cws(): s2 must contain at least one value.";
      return(.);
   end;

   if (nrow(nu) ^= n) | (nrow(w) ^= n) then do;
      print "cws(): s2, nu, and w must have the same length.";
      return(.);
   end;

   if (sum(miss(s2)) > 0) |
      (sum(miss(nu)) > 0) |
      (sum(miss(w)) > 0) then do;
      print "cws(): missing values are not allowed.";
      return(.);
   end;

   if min(s2) < 0 then do;
      print "cws(): all s2 values must be non-negative.";
      return(.);
   end;

   if min(nu) <= 0 then do;
      print "cws(): all nu values must be positive.";
      return(.);
   end;

   if min(w) <= 0 then do;
      print "cws(): all w values must be positive.";
      return(.);
   end;

   num = sum(w # s2);
   den = sum((w # w) # (s2 # s2) #/ (nu + 2));

   if den <= 0 then return(.);
   return(num**2 / den - 2);
finish cws;


start cws2(s2_1, s2_2, nu_1, nu_2, w1, w2);
   v1 = colvec(s2_1);
   v2 = colvec(s2_2);
   u1 = colvec(nu_1);
   u2 = colvec(nu_2);
   a1 = colvec(w1);
   a2 = colvec(w2);

   if (nrow(v1) ^= 1) | (nrow(v2) ^= 1) |
      (nrow(u1) ^= 1) | (nrow(u2) ^= 1) |
      (nrow(a1) ^= 1) | (nrow(a2) ^= 1) then do;
      print "cws2(): all arguments must be scalars.";
      return(.);
   end;

   if (sum(miss(v1)) > 0) | (sum(miss(v2)) > 0) |
      (sum(miss(u1)) > 0) | (sum(miss(u2)) > 0) |
      (sum(miss(a1)) > 0) | (sum(miss(a2)) > 0) then do;
      print "cws2(): missing values are not allowed.";
      return(.);
   end;

   s2_1 = sum(v1);
   s2_2 = sum(v2);
   nu_1 = sum(u1);
   nu_2 = sum(u2);
   w1   = sum(a1);
   w2   = sum(a2);

   if (s2_1 < 0) | (s2_2 < 0) then do;
      print "cws2(): s2_1 and s2_2 must be non-negative.";
      return(.);
   end;

   if (nu_1 <= 0) | (nu_2 <= 0) then do;
      print "cws2(): nu_1 and nu_2 must be positive.";
      return(.);
   end;

   if (w1 <= 0) | (w2 <= 0) then do;
      print "cws2(): w1 and w2 must be positive.";
      return(.);
   end;

   s2 = w1 * s2_1 + w2 * s2_2;
   den = w1 * w1 * s2_1 * s2_1 / (nu_1 + 2) +
         w2 * w2 * s2_2 * s2_2 / (nu_2 + 2);

   if den <= 0 then return(.);
   return(s2 * s2 / den - 2);
finish cws2;

////
//   Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
//
//   File:      djc1.java
//
//   Purpose: Demonstrates how to solve the problem with two disjunctions:
//
//      minimize    2x0 + x1 + 3x2 + x3
//      subject to   x0 + x1 + x2 + x3 >= -10
//                  (x0-2x1<=-1 and x2=x3=0) or (x2-3x3<=-2 and x1=x2=0)
//                  x0=2.5 or x1=2.5 or x2=2.5 or x3=2.5
////
package com.mosek.example;
import mosek.*;

public class djc1 {
  static double inf = 0.0; // Infinity for symbolic purposes

  public static void main (String[] args) {
    // Make an environment and task
    try (mosek.Env  env = new Env();
         mosek.Task task = new Task(env, 0, 0)) {

      // Append free variables
      int numvar = 4;
      task.appendvars(numvar);
      task.putvarboundsliceconst(0, numvar, mosek.boundkey.fr, -inf, inf);

      // The linear part: the linear constraint
      task.appendcons(1);
      task.putarow(0, new int[]{0, 1, 2, 3}, new double[]{1, 1, 1, 1});
      task.putconbound(0, mosek.boundkey.lo, -10.0, -10.0);

      // The linear part: objective
      task.putobjsense(mosek.objsense.minimize);
      task.putclist(new int[]{0, 1, 2, 3}, new double[]{2, 1, 3, 1});

      // Fill in the affine expression storage F, g
      long numafe = 10;
      task.appendafes(numafe);

      long[]   fafeidx = new long[]{0, 0, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9};
      int[]    fvaridx = new int[]{0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
      double[] fval    = new double[]{1.0, -2.0, 1.0, -3.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
      double[] g       = new double[]{1.0, 2.0, 0.0, 0.0, 0.0, 0.0, -2.5, -2.5, -2.5, -2.5};

      task.putafefentrylist(fafeidx, fvaridx, fval);
      task.putafegslice(0, numafe, g);

      // Create domains
      long zero1   = task.appendrzerodomain(1);
      long zero2   = task.appendrzerodomain(2);
      long rminus1 = task.appendrminusdomain(1);

      // Append disjunctive constraints
      long numdjc = 2;
      task.appenddjcs(numdjc);

      // First disjunctive constraint
      task.putdjc(0,                                           // DJC index
                  new long[]{rminus1, zero2, rminus1, zero2},  // Domains     (domidxlist)
                  new long[]{0, 4, 5, 1, 2, 3},                // AFE indices (afeidxlist)
                  null,                                        // Unused
                  new long[]{2, 2} );                          // Term sizes  (termsizelist)

      // Second disjunctive constraint
      task.putdjc(1,                                        // DJC index
                  new long[]{zero1, zero1, zero1, zero1},   // Domains     (domidxlist)
                  new long[]{6, 7, 8, 9},                   // AFE indices (afeidxlist)
                  null,                                     // Unused
                  new long[]{1, 1, 1, 1} );                 // Term sizes  (termidxlist)

      // Useful for debugging
      task.writedata("djc.ptf");                          // Write file in human-readable format
      // Attach a log stream printer to the task
      task.set_Stream(
        mosek.streamtype.log,
        new mosek.Stream()
      { public void stream(String msg) { System.out.print(msg); }});

      // Solve the problem
      task.optimize();

      // Print a summary containing information
      // about the solution for debugging purposes
      task.solutionsummary(mosek.streamtype.msg);

      // Get status information about the solution
      mosek.solsta solsta = task.getsolsta(mosek.soltype.itg);

      switch (solsta) {
        case integer_optimal:
          double[] xx = task.getxx(mosek.soltype.itg);

          System.out.println("Optimal primal solution\n");
          for (int j = 0; j < numvar; ++j)
            System.out.println ("x[" + j + "]:" + xx[j]);
          break;
        default:
          System.out.println("Another solution status");
          break;
      }
    }
    catch (mosek.Exception e) {
      System.out.println ("An error/warning was encountered");
      System.out.println (e.toString());
      throw e;
    }
  }
}

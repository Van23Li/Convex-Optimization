/*
   Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.

   File :      opt_server_async.java

   Purpose :   Demonstrates how to use MOSEK OptServer
               to solve optimization problem asynchronously
*/

package com.mosek.example;
import mosek.*;

public class opt_server_async {
  public static void main (String[] args) throws java.lang.Exception {
    if (args.length == 0) {
      System.out.println ("Missing argument, syntax is:");
      System.out.println ("  opt_server_async inputfile host:port numpolls");
    } else {

      String inputfile = args[0];
      String addr      = args[1];
      int numpolls     = Integer.parseInt(args[2]);
      String cert      = args.length < 4 ? null : args[3];

      try (Env env = new Env()) {
        String token;

        try(Task task = new Task(env, 0, 0)) {
          task.readdata (inputfile);
          if (cert != null)
            task.putstrparam(sparam.remote_tls_cert_path,cert);
          token = task.asyncoptimize (addr,"");
        }

        System.out.printf("Task token = %s\n", token);

        try(Task task = new Task(env, 0, 0)) {
          System.out.println("Reading input file...");

          task.readdata (inputfile);

          if (cert != null)
            task.putstrparam(sparam.remote_tls_cert_path,cert);

          System.out.println("Setting log stream...");

          task.set_Stream (mosek.streamtype.log,
          new mosek.Stream() {
            public void stream(String msg) { System.out.print(msg); }
          });

          long start = System.currentTimeMillis();

          System.out.println("Starting polling loop...");

          int i = 0;

          while ( true ) {

            Thread.sleep(100);

            System.out.printf("poll %d...\n", i);

            rescode trm[]  = new rescode[1];
            rescode resp[] = new rescode[1];

            boolean respavailable = task.asyncpoll( addr,
                                                    "",
                                                    token,
                                                    resp,
                                                    trm);


            System.out.println("polling done");

            if (respavailable) {
              System.out.println("solution available!");

              task.asyncgetresult(addr,
                                  "",
                                  token,
                                  resp,
                                  trm);

              task.solutionsummary (mosek.streamtype.log);
              break;
            }

            i++;

            if (i == numpolls) {
              System.out.println("max num polls reached, stopping host.");
              task.asyncstop (addr, "", token);
              break;
            }

          }
        } catch (java.lang.Exception e) {
          System.out.println("Something unexpected happend...");
          throw e;
        }
      }
    }
  }
}

#include "ConvexQuadOptimization.h"



bool solveConvexQuadPorgramming_mosek(std::vector<MSKboundkeye>& bkc, std::vector<double>& blc, std::vector<double>& buc,/* Bounds on constraints. */
	std::vector<MSKboundkeye>& bkx, std::vector<double>& blx, std::vector<double>& bux,/* Bounds on variables. */
	std::vector<MSKlidxt>& aptrb, std::vector<MSKidxt>& asub, std::vector<double>& aval, 
	std::vector<MSKidxt>& qsubi, std::vector<MSKidxt>& qsubj, std::vector<double>& qval, std::vector<double>& c,
	std::vector<double>& XX, bool show_success)
{
	int NUMCON = bkc.size(); int NUMVAR, NUMANZ;
	if (NUMCON == 0)
	{
		NUMVAR = c.size();
		NUMANZ = 0;
	}
	else
	{
		NUMVAR = aptrb.size() - 1;
		NUMANZ = aptrb.back();
	}
	int NUMQNZ = qval.size();

	MSKidxt       i,j;
	XX.resize(NUMVAR);

	MSKenv_t      env;
	MSKtask_t     task;
	MSKrescodee   r;

	r = MSK_makeenv(&env,NULL);
	if ( r==MSK_RES_OK )
	{
		/* Directs the log stream to the 'printstr' function. */
		MSK_linkfunctoenvstream(env,
			MSK_STREAM_LOG,
			NULL,
			printstr);
	}

	bool s = false;
	/* Initialize the environment. */   
	r = MSK_initenv(env);
	if ( r == MSK_RES_OK )
	{ 
		/* Create the optimization task. */
		r = MSK_maketask(env,NUMCON,NUMVAR,&task);

		if ( r==MSK_RES_OK )
		{
			//r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);
      
			/* Give MOSEK an estimate of the size of the input data. 
			 This is done to increase the speed of inputting data. 
			 However, it is optional. */
			if (r == MSK_RES_OK)
				r = MSK_putmaxnumvar(task,NUMVAR);
    
			if (r == MSK_RES_OK)
				r = MSK_putmaxnumcon(task,NUMCON);
      
			if (r == MSK_RES_OK)
				r = MSK_putmaxnumanz(task,NUMANZ);
  
			/* Append 'NUMCON' empty constraints.
			The constraints will initially have no bounds. */
			if ( r == MSK_RES_OK )
				r = MSK_appendcons(task, NUMCON);
  
			/* Append 'NUMVAR' variables.
			The variables will initially be fixed at zero (x=0). */
			if ( r == MSK_RES_OK )
				r = MSK_appendvars(task, NUMVAR);
  
			/* Optionally add a constant term to the objective. */
			if ( r ==MSK_RES_OK )
				r = MSK_putcfix(task,0.0);

			if (r == MSK_RES_OK)
				r = MSK_putintparam(task, MSK_IPAR_NUM_THREADS, 8);

			/*if (r == MSK_RES_OK)
			r = MSK_putdouparam(task,MSK_DPAR_INTPNT_CO_TOL_INFEAS ,1e-8);

			if (r == MSK_RES_OK)
			r = MSK_putdouparam(task,MSK_DPAR_INTPNT_CO_TOL_MU_RED ,1e-20);*/
			for(j=0; j<NUMVAR && r == MSK_RES_OK; ++j)
			{
				/* Set the linear term c_j in the objective.*/  
				if(r == MSK_RES_OK)
					r = MSK_putcj(task,j,c[j]);

				/* Set the bounds on variable j.
				blx[j] <= x_j <= bux[j] */
				if(r == MSK_RES_OK)
					r = MSK_putvarbound(task,
										j,           /* Index of variable.*/
										bkx[j],      /* Bound key.*/
										blx[j],      /* Numerical value of lower bound.*/
										bux[j]);     /* Numerical value of upper bound.*/
  
				if (NUMCON > 0)
				{
					/* Input column j of A */
					if (r == MSK_RES_OK)
						r = MSK_putacol(task,
						j,                 /* Variable (column) index.*/
						aptrb[j + 1] - aptrb[j], /* Number of non-zeros in column j.*/
						&asub[0] + aptrb[j],     /* Pointer to row indexes of column j.*/
						&aval[0] + aptrb[j]);    /* Pointer to Values of column j.*/
				}
			}
  
			/* Set the bounds on constraints.
			for i=1, ...,NUMCON : blc[i] <= constraint i <= buc[i] */
			for(i=0; i<NUMCON && r==MSK_RES_OK; ++i)
			{
				r = MSK_putconbound(task,
									i,           /* Index of constraint.*/
									bkc[i],      /* Bound key.*/
									blc[i],      /* Numerical value of lower bound.*/
									buc[i]);     /* Numerical value of upper bound.*/
			}

			if ( r==MSK_RES_OK )
			{
				/*
					* The lower triangular part of the Q
					* matrix in the objective is specified.
					*/

				/* Input the Q for the objective. */

				r = MSK_putqobj(task,NUMQNZ,&qsubi[0],&qsubj[0],&qval[0]);
			}

		
			//MSK_IPAR_INTPNT_NUM_THREADS

			if ( r==MSK_RES_OK )
			{
				MSKrescodee trmcode;

				/* Run optimizer */
				r = MSK_optimizetrm(task,&trmcode);

				/* Print a summary containing information
					about the solution for debugging purposes*/
				MSK_solutionsummary (task,MSK_STREAM_LOG);
        
				if ( r==MSK_RES_OK )
				{
					MSKsolstae solsta;
					int j;
          
					MSK_getsolsta(task,MSK_SOL_ITR,&solsta);

					switch (solsta)
					{
					case MSK_SOL_STA_OPTIMAL:
					case MSK_SOL_STA_NEAR_OPTIMAL:
						MSK_getxx(task,
							MSK_SOL_ITR,    /* Request the interior solution. */
							&XX[0]);

						if (show_success) printf("Optimal primal solution\n");
						s = true;
						break;
					case MSK_SOL_STA_DUAL_INFEAS_CER:
					case MSK_SOL_STA_PRIM_INFEAS_CER:
					case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
					case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
						printf("Primal or dual infeasibility certificate found.\n");
						break;

					case MSK_SOL_STA_UNKNOWN:
						printf("The status of the solution could not be determined.\n");
						break;
					default:
						printf("Other solution status.");
						break;
					}
				}
				else
				{
					printf("Error while optimizing.\n");
				}
			}

			if (r != MSK_RES_OK)
			{
				/* In case of an error print error code and description. */      
				char symname[MSK_MAX_STR_LEN];
				char desc[MSK_MAX_STR_LEN];

				printf("An error occurred while optimizing.\n");     
				MSK_getcodedesc (r,
									symname,
									desc);
				printf("Error %s - '%s'\n",symname,desc);
			}
		}
	}
	MSK_deletetask(&task);
	MSK_deleteenv(&env);
	return s;
}
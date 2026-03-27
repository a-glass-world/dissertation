/*************************************************************************
 *  query.c
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TRUE  1
#define FALSE 0
/*------------------------------------------------------------------------*/
int query(
          char cMsg[150], 
          int iLower, 
          int iUpper)
/*------------------------------------------------------------------------*/
	{
	int i,iHold,iVar;
	int bDone_Local;
	int neg;
	char cTemp[20];

	bDone_Local = FALSE;
	while (!bDone_Local) 
    {
		fprintf(stdout,"%s \n >> ",cMsg);
		scanf("%s",cTemp);
	   	iHold=0;
		i=0;
		iVar=0;
		neg=FALSE;
		if (cTemp[0]=='-') 
        {
			neg = TRUE;
			i++;
		}
		while ((cTemp[i]!=NULL)&&(cTemp[i]>='0')&&(cTemp[i]<='9'))	
		{
    	 	iVar=(cTemp[i]-'0')+iHold*10;   /* transform the string */
			iHold=(iVar);	 	            /* into an integer      */
			i++;
		} 
		if (neg==TRUE) iVar *= (-1);
		if (((iVar >= iLower) && (iVar <= iUpper))  /* check bounds */
		     ||(cTemp[0] == 'q'))			            /* or its ok to quit */
		{
			bDone_Local = TRUE;
			if (cTemp[0] == 'q') {
				fprintf(stderr,"...bye!\n");
				exit(0); 	
			} /* user requested exit */
		}
		else  
		{
			fprintf(stdout,"Give an integer between %d and %d.\n",iLower,iUpper);
 			fprintf(stdout,"To QUIT type \"q\"\n");
		}   /*end if-else*/
	} 
	return iVar;
} 
/*------------------------------------------------------------------------*/

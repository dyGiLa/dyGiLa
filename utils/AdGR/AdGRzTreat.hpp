#define AdGRzTreat1(A_Xmj, A_Xpj, A_Xmez, A_Xpez, abt_ratio) \
{							     \
    if ( X.coordinate(e_z) == 0 )                            \
      {                                                      \
        const real_t trcoef=((1.-(abt_ratio))/(1.+(abt_ratio))); \
	foralldir(col)                                           \
	  {                                                      \
           if(col == 2)                                          \
	     { foralldir(row){ A_Xmj.e(row, col)=0.0; } }        \
	   else                                                  \
	     { foralldir(row){ A_Xmj.e(row, col)=(A_Xpez.e(row, col))*(trcoef);} } \
	  }                                                                       \
      }                                                                           \
    else if ( X.coordinate(e_z) == (config.lz-1) )                                \
      {                                                                           \
        const real_t trcoef=((1.+(abt_ratio))/(1.-(abt_ratio)));                  \
	foralldir(col)                                                            \
	  {                                                                       \
           if(col == 2)                                                           \
	     { foralldir(row){ A_Xpj.e(row, col)=0.0; } }                         \
	   else                                                                   \
	     { foralldir(row){ A_Xpj.e(row, col)=(A_Xmez.e(row, col))*(trcoef);} } \
	  }                                                                       \
      }                                                                           \
} 

#define AdGRzTreat2(A_Xmez, A_Xpez, A_Xm_ez, A_Xp_ez, abt_ratio) \
{                                                         \
    if ( X.coordinate(e_z) == 0 )                            \
      {                                                      \
        const real_t trcoef=((1.-(abt_ratio))/(1.+(abt_ratio))); \
	foralldir(col)                                           \
	  {                                                      \
           if(col == 2)                                          \
	     { foralldir(row){ A_Xmez.e(row, col)=0.0; } }       \
	   else                                                  \
	     { foralldir(row){ A_Xmez.e(row, col)=(A_Xp_ez.e(row, col))*(trcoef);} } \
	  }                                                                       \
      }                                                                           \
    else if ( X.coordinate(e_z) == (config.lz-1) )                                \
      {                                                                           \
        const real_t trcoef=((1.+(abt_ratio))/(1.-(abt_ratio)));                  \
	foralldir(col)                                                            \
	  {                                                                       \
           if(col == 2)                                                           \
	     { foralldir(row){ A_Xpez.e(row, col)=0.0; } }                        \
	   else                                                                   \
	     { foralldir(row){ A_Xpez.e(row, col)=(A_Xm_ez.e(row, col))*(trcoef);} } \
	  }                                                                       \
      }                                                                           \
} 


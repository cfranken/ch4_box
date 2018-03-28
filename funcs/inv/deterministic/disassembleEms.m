%%% =======================================================================
%%% = disassembleEms.m
%%% = Alex Turner
%%% = 05/03/2017
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Break the emissions into a structure.
%%% =  ( 2): "assembleEms" will do the opposite.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): in -- Matrix with the emissions.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- Structure with the emissions.
%%% =======================================================================

function [ out ] = disassembleEms( in )
out.nh_ch4    = in(:,1);
out.sh_ch4    = in(:,2);
out.nh_ch4c13 = in(:,3);
out.sh_ch4c13 = in(:,4);
out.nh_mcf    = in(:,5);
out.sh_mcf    = in(:,6);
out.nh_n2o    = in(:,7);
out.sh_n2o    = in(:,8);
out.nh_c2h6   = in(:,9);
out.sh_c2h6   = in(:,10);
out.nh_oh     = in(:,11);
out.sh_oh     = in(:,12);
out.nh_co     = in(:,13);
out.sh_co     = in(:,14);
out.tau_TS    = in(:,15);
out.kX_NH     = in(:,16);
out.kX_SH     = in(:,17);
end


%%% =======================================================================
%%% = END
%%% =======================================================================

clear all
close all
clc

rule_num = 20;

  for rule = 1 : rule_num
    order_num = dunavant_order_num ( rule );
    degree = dunavant_degree ( rule );
    fprintf ( 1, '  %8d  %8d  %8d\n', rule, degree, order_num );
  end
  
  
  LEVEL_GL = 4;
  ORDER_GL = dunavant_order_num ( LEVEL_GL );
  [ z, w ] = dunavant_rule ( LEVEL_GL )
  
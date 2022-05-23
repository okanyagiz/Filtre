-- -------------------------------------------------------------
-- 
-- File Name: C:\Users\okany\Documents\MATLAB\odev_3\codegen\mlhdlc_fir\hdlsrc\mlhdlc_fir_fixpt_pkg.vhd
-- Created: 2022-05-23 23:02:58
-- 
-- Generated by MATLAB 9.11, MATLAB Coder 5.3 and HDL Coder 3.19
-- 
-- 
-- -------------------------------------------------------------


LIBRARY IEEE;
USE IEEE.std_logic_1164.ALL;
USE IEEE.numeric_std.ALL;

PACKAGE mlhdlc_fir_fixpt_pkg IS
  TYPE vector_of_signed32 IS ARRAY (NATURAL RANGE <>) OF signed(31 DOWNTO 0);
  TYPE vector_of_signed14 IS ARRAY (NATURAL RANGE <>) OF signed(13 DOWNTO 0);
  TYPE vector_of_signed37 IS ARRAY (NATURAL RANGE <>) OF signed(36 DOWNTO 0);
  TYPE vector_of_signed28 IS ARRAY (NATURAL RANGE <>) OF signed(27 DOWNTO 0);
END mlhdlc_fir_fixpt_pkg;

function [basedir,tmpfile]=gettemp(pgm)
basedir='/tmp/gettemp/';
if exist(basedir)~=7
  mkdir(basedir);
end
tmpfile=[basedir,'/',pgm];

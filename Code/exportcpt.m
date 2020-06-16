function exportcpt(cmap,range,fname,varargin)
  % EXPORTCPT(cmap,range,fname) converts the m-by-3 matrix cmap (RGB 0-1)
  % to the GMT-compatible .cpt format (RGB 0-255), with range specified by
  % range=[min,max], and writes the result to file fname. B and F values
  % are set to the min and max colors, respectively, and N is [0 0 0].
  %
  % EXPORTCPT(cmap,range,fname,B,F,N) writes the same file, but with
  % specified RGB triples (range 0-1) for B,F, and N.
  %
  % Examples:
  %   %Use default B,F,N:
  %   cmap=colormap(jet);
  %   exportcpt(cmap,[-5 5],'~/Desktop/jet.cpt');
  %
  %   %Custom B,F,N:
  %   exportcpt(cmap,[-5 5],'~/Desktop/jet.cpt',...
  %      [0 0 0],[1 1 1],[0.5 0.5 0.5])
  %
  % Eric Lindsey, Nov. 2011
  
  if(nargin == 6)
    B=varargin{1};
    F=varargin{2};
    N=varargin{3};
  elseif(nargin == 3)
      B=cmap(1,:);
      F=cmap(end,:);
      N=[0 0 0];
  else
      error(['Usage: exportcpt(cmap,range,fname) -or- '... 
        'exportcpt(cmap,range,fname,B,F,N)']);
  end
  
  m=size(cmap,1);
  
  %rescale to 0-255
  cmap=round(cmap*255);
  B=round(B*255);
  F=round(F*255);
  N=round(N*255);
  
  % open file 'fname'
  fid=fopen(fname,'w');
  
  %write header
  fprintf(fid,'%s\n','# cpt file created by: MATLAB exportcpt.m');
  fprintf(fid,'%s\n','#COLOR_MODEL = RGB');
  fprintf(fid,'%s\n','#');
  
  %write cmap values
  rng=linspace(range(1),range(2),m);
  for i=1:m-1
    fprintf(fid,'%f\t%d\t%d\t%d\t%f\t%d\t%d\t%d\n', ...
      rng(i),cmap(i,:),rng(i+1),cmap(i+1,:));
  end

  %write background, foreground, NaN specifications
  fprintf(fid,'B\t%d\t%d\t%d\n',B(1),B(2),B(3));
  fprintf(fid,'F\t%d\t%d\t%d\n',F(1),F(2),F(3));
  fprintf(fid,'N\t%d\t%d\t%d\n',N(1),N(2),N(3));
  
  fclose(fid);
  
end
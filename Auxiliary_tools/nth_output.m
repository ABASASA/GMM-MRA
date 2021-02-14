%% Extract the n-th output of thte function
function value = nth_output(N,fcn,varargin)
  [value{1:N}] = fcn(varargin{:});
  value = value{N};
end
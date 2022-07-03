function varargout = setvalues(array)
  % put array entries into several variables at once
  for i = 1:min(nargout,numel(array))
    varargout(i) = {array(i)};
  end
end

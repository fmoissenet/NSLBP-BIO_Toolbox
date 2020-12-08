classdef data_attribute < SwigRef
  methods
    function this = swig_this(self)
      this = moveckMEX(3, self);
    end
    function self = data_attribute(varargin)
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = moveckMEX(5, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.SwigClear();
      end
    end
    function delete(self)
      if self.swigPtr
        moveckMEX(6, self);
        self.SwigClear();
      end
    end
    function varargout = shape(self,varargin)
      [varargout{1:nargout}] = moveckMEX(7, self, varargin{:});
    end
    function varargout = read(self,varargin)
      [varargout{1:nargout}] = moveckMEX(8, self, varargin{:});
    end
    function varargout = write(self,varargin)
      [varargout{1:nargout}] = moveckMEX(9, self, varargin{:});
    end
  end
  methods(Static)
  end
end

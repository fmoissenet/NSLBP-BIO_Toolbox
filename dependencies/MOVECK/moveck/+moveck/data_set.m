classdef data_set < SwigRef
  methods
    function this = swig_this(self)
      this = moveckMEX(3, self);
    end
    function self = data_set(varargin)
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = moveckMEX(10, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.SwigClear();
      end
    end
    function delete(self)
      if self.swigPtr
        moveckMEX(11, self);
        self.SwigClear();
      end
    end
    function varargout = name(self,varargin)
      [varargout{1:nargout}] = moveckMEX(12, self, varargin{:});
    end
    function varargout = list_attributes_name(self,varargin)
      [varargout{1:nargout}] = moveckMEX(13, self, varargin{:});
    end
    function varargout = create_attribute(self,varargin)
      [varargout{1:nargout}] = moveckMEX(14, self, varargin{:});
    end
    function varargout = retrieve_attribute(self,varargin)
      [varargout{1:nargout}] = moveckMEX(15, self, varargin{:});
    end
    function varargout = shape(self,varargin)
      [varargout{1:nargout}] = moveckMEX(16, self, varargin{:});
    end
    function varargout = read(self,varargin)
      [varargout{1:nargout}] = moveckMEX(17, self, varargin{:});
    end
    function varargout = write(self,varargin)
      [varargout{1:nargout}] = moveckMEX(18, self, varargin{:});
    end
  end
  methods(Static)
  end
end

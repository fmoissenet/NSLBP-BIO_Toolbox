classdef data_group < SwigRef
  methods
    function this = swig_this(self)
      this = moveckMEX(3, self);
    end
    function self = data_group(varargin)
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = moveckMEX(19, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.SwigClear();
      end
    end
    function delete(self)
      if self.swigPtr
        moveckMEX(20, self);
        self.SwigClear();
      end
    end
    function varargout = name(self,varargin)
      [varargout{1:nargout}] = moveckMEX(21, self, varargin{:});
    end
    function varargout = list_group_children_name(self,varargin)
      [varargout{1:nargout}] = moveckMEX(22, self, varargin{:});
    end
    function varargout = list_set_children_name(self,varargin)
      [varargout{1:nargout}] = moveckMEX(23, self, varargin{:});
    end
    function varargout = list_attributes_name(self,varargin)
      [varargout{1:nargout}] = moveckMEX(24, self, varargin{:});
    end
    function varargout = exists_group(self,varargin)
      [varargout{1:nargout}] = moveckMEX(25, self, varargin{:});
    end
    function varargout = exists_set(self,varargin)
      [varargout{1:nargout}] = moveckMEX(26, self, varargin{:});
    end
    function varargout = create_attribute(self,varargin)
      [varargout{1:nargout}] = moveckMEX(27, self, varargin{:});
    end
    function varargout = retrieve_attribute(self,varargin)
      [varargout{1:nargout}] = moveckMEX(28, self, varargin{:});
    end
    function varargout = create_group(self,varargin)
      [varargout{1:nargout}] = moveckMEX(29, self, varargin{:});
    end
    function varargout = retrieve_group(self,varargin)
      [varargout{1:nargout}] = moveckMEX(30, self, varargin{:});
    end
    function varargout = create_set(self,varargin)
      [varargout{1:nargout}] = moveckMEX(31, self, varargin{:});
    end
    function varargout = retrieve_set(self,varargin)
      [varargout{1:nargout}] = moveckMEX(32, self, varargin{:});
    end
  end
  methods(Static)
  end
end

classdef data_store < SwigRef
  methods
    function this = swig_this(self)
      this = moveckMEX(3, self);
    end
    function self = data_store(varargin)
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = moveckMEX(33, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.SwigClear();
      end
    end
    function delete(self)
      if self.swigPtr
        moveckMEX(34, self);
        self.SwigClear();
      end
    end
    function varargout = name(self,varargin)
      [varargout{1:nargout}] = moveckMEX(35, self, varargin{:});
    end
    function varargout = list_group_children_name(self,varargin)
      [varargout{1:nargout}] = moveckMEX(36, self, varargin{:});
    end
    function varargout = list_set_children_name(self,varargin)
      [varargout{1:nargout}] = moveckMEX(37, self, varargin{:});
    end
    function varargout = list_attributes_name(self,varargin)
      [varargout{1:nargout}] = moveckMEX(38, self, varargin{:});
    end
    function varargout = exists_group(self,varargin)
      [varargout{1:nargout}] = moveckMEX(39, self, varargin{:});
    end
    function varargout = exists_set(self,varargin)
      [varargout{1:nargout}] = moveckMEX(40, self, varargin{:});
    end
    function varargout = dump(self,varargin)
      [varargout{1:nargout}] = moveckMEX(41, self, varargin{:});
    end
    function varargout = create_attribute(self,varargin)
      [varargout{1:nargout}] = moveckMEX(42, self, varargin{:});
    end
    function varargout = retrieve_attribute(self,varargin)
      [varargout{1:nargout}] = moveckMEX(43, self, varargin{:});
    end
    function varargout = create_group(self,varargin)
      [varargout{1:nargout}] = moveckMEX(44, self, varargin{:});
    end
    function varargout = retrieve_group(self,varargin)
      [varargout{1:nargout}] = moveckMEX(45, self, varargin{:});
    end
    function varargout = create_set(self,varargin)
      [varargout{1:nargout}] = moveckMEX(46, self, varargin{:});
    end
    function varargout = retrieve_set(self,varargin)
      [varargout{1:nargout}] = moveckMEX(47, self, varargin{:});
    end
  end
  methods(Static)
  end
end

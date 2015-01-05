%% Scripts used to check if there is any error in parameter configuration.
if isempty(GMM_birth)
   disp('!!WARNING: There is no birth target!'); 
else
   Dim = [];
   for n=1:length(GMM_birth)       
       nMean = size(GMM_birth(n).mean,1);
       nVar  = size(GMM_birth(n).variance,1);
       Dim = [Dim nMean]; %#ok<*AGROW>
       if (nMean~=nVar)
           error(['??ERROR:GMM_birth(' num2str(n) ').mean doesn''t match '...
                  'the dimension with GMM_birth(' num2str(n) ').variance.']);
       end
   end
   dim0 = min(Dim);
   dim1 = max(Dim);
   if dim0~=dim1
       Dimbin = dim0:dim1;
       Counter = histc(Dim, Dimbin);
       [tmp,idx1] = max(Counter); %#ok<*ASGLU>
       [tmp,idx0] = min(Counter);
       dim4 = Dimbin(idx1);
       dim3 = Dimbin(idx0);
       Idx1 = find(Dim==dim4);
       Idx0 = find(Dim==dim3,1);
       str = [];
       for m=1:length(Idx1)
           str = [str num2str(Idx1(m)) ','];
       end
       str(end) = [];
       error(['??ERROR:GMM_birth(' num2str(Idx0) ').mean doesn''t match '...
                  'the dimension with GMM_birth(' str ').mean.']);
   end
end

if ~exist('GMM_spawn','var') 
   disp('!!WARNING: There is no spawning target!'); 
else
    if isempty(GMM_spawn)
        disp('!!WARNING: There is no spawning target!');
    else
        
        Dim = [];
        for n=1:length(GMM_spawn)
            nMean = size(GMM_spawn(n).mean,1);
            nVar  = size(GMM_spawn(n).variance,1);
            Dim = [Dim nMean]; %#ok<*AGROW>
            if (nMean~=nVar)
                error(['??ERROR:GMM_spawn(' num2str(n) ').mean doesn''t match '...
                    'the dimension with GMM_spawn(' num2str(n) ').variance.']);
            end
        end
        dim0 = min(Dim);
        dim1 = max(Dim);
        if dim0~=dim1
            Dimbin = dim0:dim1;
            Counter = histc(Dim, Dimbin);
            [tmp,idx1] = max(Counter); %#ok<*ASGLU>
            [tmp,idx0] = min(Counter);
            dim4 = Dimbin(idx1);
            dim3 = Dimbin(idx0);
            Idx1 = find(Dim==dim4);
            Idx0 = find(Dim==dim3,1);
            str = [];
            for m=1:length(Idx1)
                str = [str num2str(Idx1(m)) ','];
            end
            str(end) = [];
            error(['??ERROR:GMM_spawn(' num2str(Idx0) ').mean doesn''t match '...
                'the dimension with GMM_spawn(' str ').mean.']);
        end
   end
end

if isempty(Track)
    error('??ERROR: There is no Track struct initilized!');
else
    for n=1:length(Track)
        if (Track(n).t0(1)==1) && (Track(n).isactive(1)==0)
            error(['??ERROR:Track(' num2str(n) ').t0==1 doesn''t match '...
                  'that Track(' num2str(n) ').isactive=0.']);
        end
        if (Track(n).t0(1)>1) && (Track(n).isactive(1)==1)
            error(['??ERROR:Track(' num2str(n) ').t0>1 doesn''t match '...
                  'that Track(' num2str(n) ').isactive=1.']);
        end
        if (Track(n).tf(1)>nSims) 
            error(['??ERROR:Track(' num2str(n) ').tf>' num2str(nSims) ...
                  ' which is the maximum allowable time steps.']);
        end
    end
end
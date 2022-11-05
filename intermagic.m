function varargout = intermagic(varargin)
% INTERMAGIC load "Intermagnet" magnetic observatory definitive data
%   INTERMAGIC loads time series data from the CD files of International
%       Real-time Magnetic Observatory Network.
%
% intermagic('station','series',dates);
%   where 'series' can be any of the values, 'day', 'minute', or 'hour' and
%   'station' is the three letter station code and dates is a matrix
%   containing matlab-compatable dates. (in string form, matlab date form,
%   or matlab date-vector form.)
%
%   Use this form to load from CDs as it will prompt for additional
%   directories until specified.
%
% intermagic('station','series',dates,'rootdir');
%   rootdir = path to main tree containing data files.  e.g. 'data' where
%   the relative (or full) path to a datafile is e.g.
%   'data/mag2000/abg/abg00jan.bin'
%
% intermagic('station','series',dates,'rootdir','errordir');
%   if additional files must be loaded for error corrections, specify the
%   location of those files here.  Unless specified this directory is
%   assumed to be 'rootdir/errata/'  Currently only used for loading
%   corrections to 91 data.
%
% DATA = intermagic('station','series',dates,'rootdir');
%   if output variable is specified, data will be loaded into that instead
%   of placed in the workspace in a variable named for the station ID.
%
% [DATA coords] = intermagic('station','series',dates,'rootdir');
%   if two output variables are specified, the second contains the
%   coordinate system of the stored data.  All data is returned in XYZ
%   coordinates.  If stored in HDZ, data is converted before being
%   returned.
%
% If no variables are passed as inputs, the intermagic returns the list
% of all stations with data available.
%
% example usage:
%
% DATA = intermagic('ABG','day',['8-dec-1999'; '5-mar-2000'],'data');
% DATA = intermagic('BOU','minute',['8-dec-1999'; '5-mar-2000'],'data');
%   example code to test the function.
% F = sqrt(DATA(:,1).^2 + DATA(:,2).^2 + DATA(:,3).^2); plot(DATA(:,5),DATA(:,4),DATA(:,5),F), datetick
%
% Updated 1/2006, James Wanliss

cdrom = logical(0); % If it runs from a DVD

switch nargin
    case 5
        station = varargin{1};
        series = varargin{2};
        dates = datenum(varargin{3});
        directory = varargin{4};
        erratadir =  varargin{5};
    case 4
        station = varargin{1};
        series = varargin{2};
        dates = datenum(varargin{3});
        directory = varargin{4};
        erratadir = [directory '/ERRATA/'];
    case 3
        station = varargin{1};
        series = varargin{2};
        dates = datenum(varargin{3});
        directory = 'data';
        cdrom = logical(1);
    otherwise
        printstationlist;
        return;
end

if dates(1) > dates(2), warning('Dates not consecutive.  Switching'), dates(1:2) = dates([2 1]);end

dv = datevec(dates);
coords=cellstr('NULL');

DATA = pre(diff(dates),series);
DATA(:,5) = DATA(:,5)+dates(1);

done = logical(0);
while not(done)
    if cdrom;
        dir_old = directory;
        directory = input('Enter top-level directory containing files [data/] (d=done): ','s');
        if isempty(directory)
            directory = dir_old;
        elseif strcmp(lower(directory),'d')
            break;
        end
        erratadir = [directory '/ERRATA/'];
    else, done = 1;
    end
    filelist = localfind([directory '/'],['*.bin'],[station '/']);
    for j = 1:length(filelist)
        [start finish]  = regexpi(filelist(j).name,'\d{4}');
        fileyear = str2num(filelist(j).name([start:finish]))
        if (fileyear <= dv(2,1)) && (fileyear >= dv(1,1))
            [start finish]  = regexpi(filelist(j).name,station);
            try
                if(strcmp(station,filelist(j).name(start(1):finish(1))))
                    [DATA_M coords]= loadfile(filelist(j),series);
                    DATA = wisemerge(DATA,DATA_M,filelist(j).name);
                    PlotCompleted(DATA(:,1:4),DATA(:,5))
                end
            catch
                disp(strcat('Do not load ',filelist(j).name));
            end
        end
    end
end

DATA = fixdata(DATA,station,erratadir,series);
disp(['Orientation: ' coords{1}])
if strfind(coords{1},'D'),
    DATA(:,1:3) = HDZtoXYZ(DATA(:,1:3),coords{1});
end

switch nargout
    case 0
        assignin('base',upper(station),DATA)
    case 1
        varargout = {DATA};
    case 2
        varargout = {DATA ;coords};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotCompleted(S,T);
% adapted from matlab "spy" function.  markers indicate where data is still
% missing.
[m,n] = size(S);
[i,j] = find(isnan(S));
if isempty(i), i = NaN; j = NaN; else i = T(i); end

plot(i,j,'marker','o','linestyle','none','markersize',4,'color','r')

[i,j] = find(~isnan(S));
if isempty(i), i = NaN; j = NaN; else i = T(i); end
hold on
plot(i,j,'marker','.','markersize',2,'linestyle','none')
hold off
if ~isnan(i), datetick('x'); end

axis tight
Title('Circled elements remain to be filled')
set(gca,'ylim',[0 n+1],'ydir','reverse', ...
    'grid','none','plotboxaspectratio',[10 1 1],...
    'ytick',1:4,'yticklabel',{'X' 'Y' 'Z' 'F'});


function DATA = fixdata(DATA,station,erratadir,series)
DATA(DATA==999999) = NaN;
DATA(:,1:4) = DATA(:,1:4)*.1; %convert 10s of nT to nT.
DV = datevec([DATA(1,5);DATA(size(DATA,1),5)]);
STADT = strcat(upper({station}), num2str(mod(DV(1,1):DV(2,1),100)','%0.2d'));


for i = 1:length(STADT)
    switch STADT{i}
        % Known errors in INTERMAGNET CD-ROM's (from 2003 CD)
        % When updates are available, paste the changes as comments at the
        % end of this file and add code making the changes to an
        % appropriate case.
        %%%%%%%%%%%%%%%%%%%%%%%%% 1991 CD-ROM:
        % Errors in the hourly mean values for the Y-components for the following
        % observatories: BLC, CBB, FCC, MBC, MEA, OTT, STJ and VIC. The correct values
        % are given in files on 1994 CD-ROM.
        case {'BLC91', 'CBB91', 'FCC91', 'MBC91', 'MEA91', 'OTT91', 'STJ91', 'VIC91'}
            while ~((exist([erratadir STADT{i} 'Y.MHV']) == 2) || strcmp(erratadir,'s'))
                erratadir = input(['Enter directory containing ''' STADT{i} 'Y.MHV'' (s=skip): '],'s');
            end
            if strcmp(erratadir,'s'),
                warning('ERRATA:uncorrected','1991: hourly mean values for the Y-components for station: %s not corrected',STADT{i}(1:3));
                continue;
            end;
            warning('ERRATA:corrected','1991: Errors corrected in the hourly mean values \nfor the Y-components for station: %s',STADT{i}(1:3))
            RepData = load91err([erratadir STADT{i} 'Y.MHV']);
            [c,ia,ib] = intersect(DATA(:,5),RepData(2,:));
            DATA(ia,2) = RepData(1,ib)';


        case {'CLF91', 'NUR91', 'SOD91'}
            % Ak numbers in binary files are multiplied by 10 for CLF, NUR, SOD.
            warning('ERRATA:uncorrected','Not Corrected: 1991:  Ak numbers in binary files are multiplied by 10 for CLF, NUR, SOD.')

        case {'PPT91', 'THY91'}
            % Changes with replacement files:
            % 	PPT: Errors in recorded hourly means. Scaling error in D-component. Replace-
            %          ment files on 1998 CD-ROM.
            % 	THY: Errors in daily means on day 37, 38, 235 and 267. Replacement files on
            %          1998 CD-ROM.
            warning('ERRATA:corrected','Corrected: 1991: Station %s errors fixed in replacement file',station)
        case 'TIK91'
            % 	TIK: Recorded hourly means for Z-component in April are 10 times too small.
            warning('ERRATA:corrected','TIK: Recorded hourly means for Z-component in April were 10 times too small.')
            if strcmp(series,'hour')
                y = datenum(1991,4,0,0,0,0); z = datenum(1991,5,0,0,0,0);
                DATA((DATA(:,5)>=y) & (DATA(:,5)<z),3 ) = DATA((DATA(:,5)>=y) & (DATA(:,5)<z),3 )*10;
                clear y z
            end
        case 'BSL91'
            % BSL: Large Data spike in H on June 4.
            warning('ERRATA:nocorrection','1991 station BSL: Large Data spike in H on June 4.')
            %
            %%%%%%%%%%%%%%%%%%%%%%%%% 1992 CD-ROM:

        case {'NUR92', 'SOD92', 'DRV92', 'THY92'}
            % Ak numbers in binary files are multiplied by 10 for CLF, NUR, SOD, DRV, THY.
            warning('ERRATA:uncorrected','Ak Numbers in binary files incorrect')
        case {'CLF92'}
            % CLF: Hourly means are not correct.
            warning('ERRATA:uncorrected','Ak Numbers in binary files incorrect')
            warning('ERRATA:corrected','Hourly means are not correct, replaced with NaN')
            if strcmp(series,'hour')
                y = datenum(1992,0,0,0,0,0); z = datenum(1993,0,0,0,0,0);
                DATA((DATA(:,5)>=y) & (DATA(:,5)<z),1:4) = NaN;
                clear y z
            end
            %
            %%%%%%%%%%%%%%%%%%%%%%%%% 1993 CD-ROM:
        case {'CLF93', 'NUR93', 'SOD93', 'DRV93', 'THY93'}
            % Ak numbers in binary files are multiplied by 10 for CLF, NUR, SOD, DRV, THY.
            % Ak numbers in binary files for NCK are all 0.
            warning('ERRATA:uncorrected','Ak numbers uncorrected, but unused')
        case {'NCK93'}
            % NCK: D-baseline values are given in 1/10 nT not in 1/10 minute of arc. The
            % scale of the baseline plot for D is therefore a factor of 6 too big.
            warning('ERRATA:uncorrected','Correction not made to D-baseline value as it is unused')
        case {'FRN93','HON93'}
            % FRN: Errors in K-numbers in binary files. Replacement data on 2003 CD-ROM.
            % HON: Errors in K-numbers in binary files. Replacement data on 2003 CD-ROM.
            warning('ERRATA:uncorrected','Correction not made to K-numbers')
            %
            %%%%%%%%%%%%%%%%%%%%%%%%% 1994 CD-ROM: %%%%%%%%%%%%%%%%%%%%%%%%
        case {'CLF94', 'DRV94', 'THY94'}
            % Ak numbers in binary files are multiplied by 10 for CLF, NUR, SOD, DRV, THY.
            % Ak numbers in binary files for AMS are all 0.
            warning('ERRATA:uncorrected','Correction not made to ak numbers')
        case {'EYR94'}
            % EYR: Longitude and Co-Latitude must be multiplied by 100
            warning('ERRATA:uncorrected','Correction not made to longitude & co-latitude')
        case {'HER94','WNG94'}
            % HER: Country map missing
            % WNG: Missing * at the end of BLV-file
            warning('ERRATA:nondata','Errors in non-data files')
        case {'NUR94','SOD94'}
            % NUR: F-component is wrong, should not be used.
            % SOD: F-component is wrong, should not be used.
            warning('ERRATA:uncorrected','Correction not made to ak numbers')
            warning('ERRATA:corrected','F-component wrong - changed to NaN')
            y = datenum(1994,0,0,0,0,0); z = datenum(1995,0,0,0,0,0);
            DATA((DATA(:,5)>=y) & (DATA(:,5)<z),4) = NaN;
            clear y z

            %
            %%%%%%%%%%%%%%%%%%%%%%%%% 1995 CD-ROM: %%%%%%%%%%%%%%%%%%%%%%%%

        case {'CLF95','NUR95', 'SOD95', 'DRV95', 'THY95'}
            % Ak numbers in binary files are multiplied by 10 for CLF, NUR, SOD, DRV, THY.
            % Ak numbers in binary files for NCK are not correct.
            warning('ERRATA:uncorrected','Ak numbers incorrect')

        case {'EYR95'}
            % EYR: Longitude and Co-Latitude must be multiplied by 100
            warning('ERRATA:uncorrected','Not  Corrected: 1995 Longitude and Co-Latitude')
        case {'NUR95','SOD95'}
            % NUR: F-component is wrong, should not be used.
            % SOD: F-component is wrong, should not be used.
            warning('ERRATA:corrected','F-component wrong - changed to NaN')
            y = datenum(1995,0,0,0,0,0); z = datenum(1996,0,0,0,0,0);
            DATA((DATA(:,5)>=y) & (DATA(:,5)<z),4) = NaN;
            clear y z

        case {'NCK95'}
            % NCK: D-baseline values are given in 1/10 nT not in 1/10 minute of arc. The
            % scale of the baseline plot for D is therefore a factor of 6 too big.
            %
            warning('ERRATA:uncorrected','Not Corrected: 1995: D-baseline values in 1/10 nT')
            %
            %%%%%%%%%%%%%%%%%%%%%%%%% 1996 CD-ROM: %%%%%%%%%%%%%%%%%%%%%%%%

        case {'CLF96', 'NUR96', 'SOD96', 'DRV96'}
            % Ak numbers in binary files are multiplied by 10 for CLF, NUR, SOD, DRV, THY.

            warning('ERRATA:uncorrected','Not Corrected: 1996: ak numbers incorrect')
        case {'EYR96', 'SBA96'}
            % EYR and SBA: Longitude and Co-Latitude must be multiplied by 100
            %             warning('ERRATA:uncorrected','Not yet fixed: 1996: EYR and SBA: Longitude and Co-Latitude must be multiplied by 100')

            warning('ERRATA:uncorrected','Correction not made: longitude and co-latitude')
        case {'THY96'}
            % THY: year in binary files is 1997!! %corrected.
            warning('ERRATA:uncorrected','Not Yet Corrected: 1996: ak numbers incorrect')

        case {'CLF96'}
            % CLF: D-component is wrong. Replacement files on 1998 CD-ROM.
            warning('ERRATA:corrected','Corrected: 1996: D-component.  Corrections in replacement files')
        case {'NUR96','SOD96'}
            % NUR: F-component is wrong, should not be used.
            % SOD: F-component is wrong, should not be used.
            warning('ERRATA:corrected','F-component wrong - changed to NaN')
            y = datenum(1996,0,0,0,0,0); z = datenum(1997,0,0,0,0,0);
            DATA((DATA(:,5)>=y) & (DATA(:,5)<z),4) = NaN;
            clear y z

        case {'NCK96'}
            % NCK: D-baseline values are given in 1/10 nT not in 1/10 minute of arc. The
            % scale of the baseline plot for D is therefore a factor of 6 too big.
            warning('ERRATA:uncorrected','Correction not made D-baseline value')
            %
            %%%%%%%%%%%%%%%%%%%%%%%%% 1997 CD-ROM: %%%%%%%%%%%%%%%%%%%%%%%%

        case {'NUR97','SOD97'}
            % NUR: F-component is wrong, should not be used.
            % SOD: F-component is wrong, should not be used.
            warning('ERRATA:corrected','F-component wrong - changed to NaN')
            y = datenum(1997,0,0,0,0,0); z = datenum(1998,0,0,0,0,0);
            DATA((DATA(:,5)>=y) & (DATA(:,5)<z),4) = NaN;
            clear y z
        case {'NCK97'}
            % NCK: D-baseline values are given in 1/10 nT not in 1/10 minute of arc. The
            % scale of the baseline plot for D is therefore a factor of 6 too big.
            warning('ERRATA:uncorrected','Correction not made: D-baseline')

            %
            %%%%%%%%%%%%%%%%%%%%%%%%% 1998 CD-ROM: %%%%%%%%%%%%%%%%%%%%%%%%
        case {  'AAE98'    'ABG98'    'ABK98'    'AMS98'    'BDV98'    'BEL98'    'BFE98'    'BLC98',...
                'BMT98'    'BNG98'    'BOU98'    'BRW98'    'BSL98'    'CBB98'    'CLF98'    'CMO98',...
                'CNB98'    'CZT98'    'DLR98'    'DRV98'    'ESK98'    'EYR98'    'FCC98'    'FRD98',...
                'FRN98'    'FUR98'    'GDH98'    'GNA98'    'GUA98'    'GUI98'    'HAD98'    'HER98',...
                'HLP98'    'HON98'    'HRB98'    'IQA98'    'IRT98'    'KAK98'    'KOU98'    'LER98',...
                'LNP98'    'LOV98'    'MBO98'    'MEA98'    'MMB98'    'NAQ98'    'NCK98'    'NEW98',...
                'NGK98'    'NUR98'    'OTT98'    'PAF98'    'PBQ98'    'PHU98'    'PPT98'    'RES98',...
                'SBA98'    'SIT98'    'SJG98'    'SOD98'    'SPT98'    'STJ98'    'TAM98'    'TAN98',...
                'THL98'    'THY98'    'TUC98'    'VIC98'    'WNG98'    'YKC98'}
            % Ak numbers on 1998 CD-ROM should not be used, since these numbers may not
            % be correct. Ak numbers are no longer published on the CD-ROM's.
            warning('ERRATA:uncorrected','Ak numbers on 1998 CD-ROM should not be used\n\tsince these numbers may not be correct.\nAk numbers are no longer published on the CD-ROM''s.\n')
        case {'NCK98'}
            % NCK: D-baseline values are given in 1/10 nT not in 1/10 minute of arc. The
            % scale of the baseline plot for D is therefore a factor of 6 too big.
            warning('ERRATA:uncorrected','1998: D baseline values in 1/10 nT')
            %
            %%%%%%%%%%%%%%%%%%%%%%%%% 1999 CD-ROM %%%%%%%%%%%%%%%%%%%%%%%%
        case {'VSS99'}
            % VSS: F-component shifted 200 nT in June
            warning('ERRATA:uncorrected','not yet corrected: 1999: F-component shifted 200nT in June')
        case {'NCK99'}
            % NCK: D-baseline values are given in 1/10 nT not in 1/10 minute of arc. The
            % scale of the baseline plot for D is therefore a factor of 6 too big.
            warning('ERRATA:uncorrected','not yet corrected: 1999: D-baseline values in nT')
            %
            %%%%%%%%%%%%%%%%%%%%%%%%% 2000 CD-ROM %%%%%%%%%%%%%%%%%%%%%%%%

        case {'ABG00'}
            % ABG: Errors in daynumbers MAR-DEC 2000. Replacement files on 2001 CD-ROM.
            % South Africa: About screen in wrong format!
            warning('ERRATA:corrected','Corrected: 2000: Station %s errors fixed in replacement file',STADT{i}(1:3))
        case {'VSS00'}
            % VSS: Data for Oct.8 and Oct.9 are identical. Data for Oct.8 are wrong and
            % should not be used.
            warning('ERRATA:uncorrected','not yet corrected: 2000: oct 8 data incorect')
            y = datenum(2000,10,8,0,0,0); z = datenum(2000,10,9,0,0,0);
            DATA((DATA(:,5)>=y) & (DATA(:,5)<z),1:4) = NaN;
            clear y z


            %%%%%%%%%%%%%%%%%%%%%%%%% 2001 CD-ROM %%%%%%%%%%%%%%%%%%%%%%%%

        case {'TAN01','HONO1'}
            % TAN: X- and Y data not correct. Also errors in yearmean-file.
            % Replacement data on 2002 CD-ROM.
            % HON:Minor errors in baseline computation
            % Replacement data on 2003 CD-ROM.
            warning('ERRATA:corrected','Corrected: 2001: Station %s errors fixed in replacement file',STADT{i}(1:3))
            %
            %%%%%%%%%%%%%%%%%%%%%%%%% 2002 CD-ROM %%%%%%%%%%%%%%%%%%%%%%%%

        case {'THY02'}
            % THY: K-numbers in binary files are not correct.
            warning('ERRATA:nocorrection','No correction: 2002: THY K-numbers incorrect')
        case {'HON02','PPT02'}
            % HON:Minor errors in baseline computation
            % Replacement data on 2003 CD-ROM.
            % PPT: incorrect jump on 2002. Error on Azimuth marks of new pillar.
            % Replacement data on 2003 CD-ROM.
            warning('ERRATA:corrected','Corrected: 2001: Station %s errors fixed in replacement file',STADT{i}(1:3))
    end
end

function DATA = HDZtoXYZ(DATA,coords)
warning('INTERMAGIC:orientation','Orientation: %s changed to XYZF',coords)
D = strfind(coords,'D'); H = strfind(coords,'H');
Z = strfind(coords,'Z');

declination = DATA(:,D)/60*pi/180; %radians
% X = Hcos(D), Y = Hsin(D)
X = DATA(:,H).*cos(declination); Y = DATA(:,H).*sin(declination);

DATA(:,1:3) = [X, Y, DATA(:,Z) ];

function DATA = wisemerge(DATA,DATA_M,f)
% merge data from the monthly data array into the main data array.  If the
% data already in the main array was put their by a "correction" file, do
% not replace it.  If the monthly data is from a "correction" file, replace
% data in the main array.

[c,ia,ib] = intersect(DATA(:,5),DATA_M(:,5));

if isempty(ia)
    %do nothing
elseif min(min(isnan(DATA(ia,1:4))))
    DATA(ia,1:4) = DATA_M(ib,1:4);
elseif isempty(strfind(upper(f),'MAG'))
    DATA(ia,1:4) = DATA_M(ib,1:4);
end


function files = localfind(directory,match,subdir)
% recursively find all files below "directory" matching the match terms.
% for speed, skips all directories of 3 letters or less (excludes the
% majority of station directories except the one of interest: subdir and
% the "up one dir" and "this dir" identifiers: '.' and '..' )
disp(['searching directory: ',directory]);
contents = dir(directory);
files = dir([directory subdir match]);
for j = 1:length(files)
    files(j).name = [subdir files(j).name];
end

files = [files;  dir([directory match])];

for j = 1:length(files)
    files(j).name = [directory files(j).name];
end

for j = 1:length(contents)
    if (contents(j).isdir) && (length(contents(j).name) > 3) && isempty(regexpi(directory,'mag'))
        ff = localfind([directory contents(j).name '/'],match,subdir);
        files = [files; ff];
    end
end


function [DATA, orientation] = loadfile(f,series)
% read series data from file specified.  return the data and the coordinate
% system used.
disp(['loading from file: ',f.name])
switch series
    case 'day'
        fieldsize = 1;
        startpos = 5872;
    case 'hour'
        fieldsize = 24;
        startpos = 5776;
    case 'minute'
        fieldsize = 1440;
        startpos = 16;
end
RecordsToRead = f.bytes/(5888*4);
DATA = zeros(RecordsToRead*fieldsize,5);

fid = fopen(f.name,'r');

% Check that orientation is homogenous.
fseek(fid,5*4,-1);
orientation = unique(cellstr(fread(fid,[4,RecordsToRead],[ '4*char=>char'],(5888-1)*4)'));
if length(orientation) > 1, error('INTERMAGIC:loaddata','Orientation not homonogeous: %s,%s',orientation{:}), end

% Check that station ID is homogenous.
fseek(fid,0,-1);
station = unique(cellstr(fread(fid,[4,RecordsToRead],[ '4*char=>char'],(5888-1)*4)'));
if length(station) > 1, error('INTERMAGIC:loaddata','stationID not homonogeous: %s,%s',station{:}), end


% read file once, parse data into array
fseek(fid,startpos*4,-1);
X = fread(fid,[4*fieldsize,RecordsToRead],[int2str(fieldsize*4) '*int32'],(5888-fieldsize*4)*4);
for component = 0:3
    Y = X((1:fieldsize)+fieldsize*component,1:RecordsToRead);
    DATA(:,component+1) = Y(:);
end
fseek(fid,1*4,-1);
T = fread(fid,RecordsToRead,'1*int32',(5888-1)*4);

fclose(fid);

Ts = datevec(datenum(floor(T/1000),0,mod(T,1000),0,0,0));

% THY96: year in binary files is 1997!!
if strfind(f.name,'THY96'),
    Ts(:,1) = 1996;
    warning('ERRATA:corrected','THY year changed back to 1996')
end
Ts = datenum(Ts);
[a b] = meshgrid(Ts,(0:fieldsize-1)/fieldsize);
z = a+b;
DATA(:,5) = z(:);

function EM = pre(sz,series)
%preallocate data array
switch series
    case 'day'
        fieldsize = 1;
    case 'hour'
        fieldsize = 24;
    case 'minute'
        fieldsize = 1440;
end

EM = repmat(NaN,sz*fieldsize+1,5);
EM(:,5)=(0:sz*fieldsize)'/fieldsize;

function printstationlist
disp('Stations available: ')
disp({...
    'AAE'    'CMO'    'HAD'    'MCQ'    'SOD'
    'ABG'    'CNB'    'HBK'    'MEA'    'SPT'
    'ABK'    'CTA'    'HER'    'MID'    'STJ'
    'ALE'    'CZT'    'HLP'    'MMB'    'SUA'
    'AMS'    'DLR'    'HON'    'NAQ'    'TAM'
    'API'    'DOU'    'HRB'    'NCK'    'TAN'
    'AQU'    'DRV'    'HRN'    'NEW'    'TEO'
    'ASC'    'EBR'    'HUA'    'NGK'    'THL'
    'ASP'    'ESK'    'IQA'    'NUR'    'THY'
    'BDV'    'EYR'    'IRT'    'OTT'    'TIK'
    'BEL'    'FCC'    'KAK'    'PAF'    'TRW'
    'BFE'    'FRD'    'KDU'    'PBQ'    'TUC'
    'BLC'    'FRN'    'KNY'    'PHU'    'UPS'
    'BMT'    'FUR'    'KOU'    'PPT'    'VAL'
    'BNG'    'GDH'    'LER'    'PST'    'VIC'
    'BOU'    'GLN'    'LNP'    'QSB'    'VSS'
    'BRW'    'GNA'    'LOV'    'RES'    'WNG'
    'BSL'    'GUA'    'LZH'    'SBA'    'YKC'
    'CBB'    'GUI'    'MBC'    'SIT'      ''
    'CLF'    'GZH'    'MBO'    'SJG'      '' })

function ou = load91err(filename)
%load y-data corrections for stations where required.
fid = fopen(filename,'r');
for i = 1:53
    (fgetl(fid));
end
S = [];
while 1
    tline = fgetl(fid);
    if isempty(tline), break, end
    S = [S;tline];
end
fclose(fid);
i = 1:size(S,1);
code = S(i,1:3); year2d = S(i,4:5); month = S(i,6:7); dayofmonth = S(i,9:10);
dates = datenum(strcat(month,'-',dayofmonth,'-',year2d));
DTA = sscanf(cell2mat(regexprep(S(i,17:120),'(.{4})',' $1','tokenize'))','%i',[26,size(S,1)]);
DTA(DTA==9999)=NaN;
F = DTA(1,:);
Xi = DTA(2:25,:);
Xi = Xi + 100*meshgrid(F,Xi(:,1));
[days hours] = meshgrid(dates,0:23);
Di = days + hours/24;
ou = [Xi(:) Di(:)]';

%%%%%%%%%%%%%%%%%%%%   Datafile specification   %%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The data on the CD-ROM are coded as 32-bit (four byte) (long integer) binary words,
%     with 5888 words comprising a day-long record. Each file contains one
%     month of day records.
%
% Words 1-10      contain header information including a
%  1:             3-letter observatory identification (ID) code
%  2:             the year concatenated with the day of the year
%                       (ex. 1999001 = jan 1, 1999)
%  3:             co-latitude
%  4:             longitude
%  5:             elevation
%  6:             orientation
%  7:             originating organization
%  8:             a D-conversion factor
%  9:             data quality
% 10:             instrumentation.
%
%                 The D-conversion factor is a fixed value used only in the
%                 graphics portion of the access software to allow
%                 Declination to be plotted in minutes of arc and
%                 equivalent nanoteslas (nT). It is given as H/3438*10000,
%                 where H is the annual mean value of the horizontal
%                 intensity. ASCII values, such as the observatory ID and
%                 orientation, are also stored as 32-bit words, but are
%                 coded as the hexadecimal byte-string corresponding to the
%                 ASCII string. For example, the string "HDZF" is coded as
%                 the sequence "48 44 5A 46".
%
% Word 11         K-9 value for observatory in nT
% Word 12         is digital sampling rate in msec
% Word 13         is sensor orientation.
%                   Sensor orientation could by XYZF, DIF, UVZ, HDZ, HDZF
%                   etc. It should indicate which components are actually
%                   measured
%
% Words 14-16     are reserved for future use and padded with zeros.
% Words 17-5776   contain the minute values of the 4 components
%                 (successively H,D,Z,F or X,Y,Z,F) for the day. The values
%                 are stored in tenth- units with an implied decimal point.
%                 Thus, an H value of 21305.6 is stored (in tenth-nT) as
%                 213056 with a decimal point implied between the last and
%                 next-to-last digits.
%
% Words 5777-5872 are used for the hourly mean values of the successive 4
%                 components.
% Words 5873-5876 store the 4 daily mean values.
% Words 5877-5884 hold the 8 K-Index values for the day.
% Words 5885-5888 set aside for each contributing institution to use as
%                 they wish, provided it is coded as a 32-bit binary value.
%
% Missing data for minute, hour, and day values are stored as "999999".
% Missing K-Index are stored as "999".


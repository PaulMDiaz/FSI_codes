function [ colors ] = get_color_palettes(testColors)

%==========================================================================
%==========================================================================
%
% Returns a structure containing lots of fun colors to use for plotting.
% Can also be used to create a test plot to see what different colors look
% like together
%
%
% Author: Eric W. Ferguson, Luke Bury
% Date: 08/26/16
%
%
% INPUT:         Description                                          Units
%  testColor     (optional)...[nx3] matrix of colors you wish to see
%                plotted
%
% OUTPUT:       
%    
%  colors        structure of colors where each field is a           struct
%                pallete and field with a pallete are RGB colors                
%
%==========================================================================
%==========================================================================
colors = struct;
%==========================================================================
%% Standard Colors
%==========================================================================
colors.std.blue       = [50  100 200]./255;
colors.std.ltblue     = [125 216 255]./255;
colors.std.orange     = [255 127 14]./255;
colors.std.ltorange   = [255 187 120]./255;
colors.std.grn        = [44 160 44]./255;
colors.std.ltgrn      = [152 223 138]./255;
colors.std.red        = [209   0   0]./255;
colors.std.ltred      = [255 152 150]./255;
colors.std.purp       = [148 103 189]./255;
colors.std.ltpurp     = [197 176 213]./255;
colors.std.brown      = [140 86 75]./255;
colors.std.ltbrown    = [196 156 148]./255;
colors.std.pink       = [227 119 194]./255;
colors.std.ltpink     = [247 182 210]./255;
colors.std.grey       = [127 127 127]./255;
colors.std.ltgrey     = [199 199 199]./255;
colors.std.ylwgrn     = [188 189 34]./255;
colors.std.khaki      = [219 219 141]./255;
colors.std.turq       = [23 190 207]./255;
colors.std.ltturq     = [158 218 229]./255;
colors.std.mag        = [230  47 215]./255;
colors.std.maglt      = [229 194 237]./255;
colors.std.ylw        = [247 202   0]./255;
colors.std.black      = [0 0 0]./255;

%==========================================================================
%% Color Schemes
%==========================================================================
%%%------------------------
%%% Risk color schemes
%%%------------------------
colors.sch.r6 = [179 19 19;240 80 83;246 144 61;247 212 32;164 224 87;25 214 29]./255;
colors.sch.r9 = [107 11 11;179 19 19;240 80 83;246 127 61;247 180 32;247 212 32;217 205 32;164 224 87;25 214 29]./255;

%%%------------------------
%%% Distinct color schemes
%%%------------------------
%%% 3 colors
colors.sch.d3_1 = [239 71 111; 0 136 204; 4 173 74]./255; % redish, blueish, greenish

%%% 4 colors
colors.sch.d4_1 = [6 214 160; 38 84 124; 239 71 111; 255 196 61]./255; % 4 distinct colors #1

%%%------------------------
%%% Similar color schemes
%%%------------------------
%%% 4 colors
colors.sch.s3_1 = [0 43 65;25 100 126;40 175 176]./255; % Blues

%==========================================================================
%% Color Testing
%==========================================================================
if nargin ~= 0
    figure
    Y = ones(2,size(testColors,1));
    h = area(Y);
    for kk = 1:size(testColors,1)
        h(kk).FaceColor = testColors(kk,:);
    end
end

end


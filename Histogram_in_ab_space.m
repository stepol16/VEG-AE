function [ Hist_ab, Hist_3D ] = Histogram_in_ab_space ( InImg, Coeff );

%
% The function evaluate a normalized (sum to 1) histogram
% of the input image in the ab plane of the CIELab color space.
% The histogram is evaluated ONLY in the region specified
% by the aLim and bLim vectors, whose values are set
% empirically.
%
% Hist_ab is a matrix whose elements are the histogram
%
% Hist_3D stores, as a structure n x 3, the elements of the histogram
% where:
%           Hist_3D ( i, 1 ) is the a coordinates of the cell
%           Hist_3D ( i, 2 ) is the b coordinates of the cell
%           Hist_3D ( i, 3 ) is the histogram value of the pair ( a, b )
%
FlagDebug = 0;
MyEps     = 1;

%
% In the TmpImg matrix n x 3 (double) the colors of
% all the non-zero pixels of the image are copied.
%
DInImg = double ( InImg );
DInImg = reshape ( DInImg, [], 3 );
if max(max( DInImg )) > 1
    DInImg = DInImg / 255.0;
end
Idx = find ( ( DInImg ( :, 1 ) ~= 0 ) | ...
             ( DInImg ( :, 2 ) ~= 0 ) | ...
             ( DInImg ( :, 3 ) ~= 0 ) );
TmpImg = DInImg ( Idx, : );

%
% TmpLab is the CIELab version of the variable TmpImg
%
TmpLab = double ( rgb2lab ( TmpImg, 'WhitePoint', 'd65' ) );

%
% The mimimum and maximum values for the two
% components a and b are evaluated
%
aMin = min ( TmpLab ( :, 2 ) );
aMax = max ( TmpLab ( :, 2 ) );
bMin = min ( TmpLab ( :, 3 ) );
bMax = max ( TmpLab ( :, 3 ) );

if FlagDebug
    fprintf ( 1, 'Min(a) = %f, Max(a) = %f, Min(b) = %f, Max(b) = %f \n', aMin, aMax, bMin, bMax );
end

%
% The limits of the area of interest in the ab plane are set
% The corresponding range is evaluated (empirically)
%
aLim = double ( [ -130.0 130.0 ] );
aRng = ceil ( double ( aLim ( 2 ) - aLim ( 1 ) ) );
bLim = double ( [ -130.0 130.0 ] );
bRng = ceil ( double ( bLim ( 2 ) - bLim ( 1 ) ) );

if ( ( aMin < aLim ( 1 ) ) | ...
     ( aMax > aLim ( 2 ) ) | ...
     ( bMin < bLim ( 1 ) ) | ...
     ( bMax > bLim ( 2 ) ) )
    fprintf ( 1, 'WARNING: there are colors out of the considered ranges \n' );
end

%
% The number of bins along each dimension of the plane
% is set (empirically)
%
% NBin = max ( aRng, bRng );
NBin = max ( aRng, bRng ) + 1;

% Reduce/group the Bins by a coefficient 
NBin = floor(NBin / Coeff);

%
% The variable X and Y are set according to the area
% of interest to enable the final plot of the histogram
% The variables a and b are associated to the axis X and Y respectively
%
X = double ( zeros ( NBin, 1 ) );
Y = double ( zeros ( NBin, 1 ) );
for i = 1 : NBin
    
    X (i) = double (i - 1 + aLim ( 1 ) );
    Y (i) = double (i - 1 + bLim ( 1 ) );

end

%
% The output variable Hist_ab is initialized
%
Hist_ab = double ( zeros ( NBin, NBin ) );
Hist_a  = double ( zeros ( NBin, NBin ) );
Hist_b  = double ( zeros ( NBin, NBin ) );
Hist_3D = double ( zeros ( NBin * NBin, 3 ) );

%
% The variable TotHist is used to normalize the histogram
% and counts all the pixels inside the area of interest
% of the ab plane
% The esclusion of the area around the zero, using the variable MyEps
% avoids the spike around zero that, due to the normalization, would reduce
% the visibility of the other components of the histogram
%
TotHist = 0;
for i = 1 : size ( TmpLab, 1 )
    
    if ( ( TmpLab ( i, 2 ) >=  MyEps ) | ...
         ( TmpLab ( i, 2 ) <= -MyEps ) | ...
         ( TmpLab ( i, 3 ) >=  MyEps ) | ...
         ( TmpLab ( i, 3 ) <= -MyEps ) )
     
         if ( ( TmpLab ( i, 2 ) >= aLim ( 1 ) ) & ...
              ( TmpLab ( i, 2 ) <= aLim ( 2 ) ) & ...
              ( TmpLab ( i, 3 ) >= bLim ( 1 ) ) & ...
              ( TmpLab ( i, 3 ) <= bLim ( 2 ) ) )

                idxa = double (  TmpLab ( i, 2 ) - aLim ( 1 ) + 1 ) / Coeff;
                idxa = int16 ( floor ( idxa ) );
                idxb = double (  TmpLab ( i, 3 ) - bLim ( 1 ) + 1 ) / Coeff;
                idxb = int16 ( floor ( idxb ) );
                Hist_ab ( idxa, idxb ) = Hist_ab ( idxa, idxb ) + 1;
                TotHist = TotHist + 1;
         end
            
    end
    
end
[ Hist_a, Hist_b ] = meshgrid ( X, Y );

%
% The variable Hist_ab is normalized (sum to 1) using
% TotHist
%
Hist_ab = Hist_ab / double ( TotHist );

Hist_3D ( :, 1 ) = reshape ( Hist_a, [], 1 );
Hist_3D ( :, 2 ) = reshape ( Hist_b, [], 1 );
Hist_3D ( :, 3 ) = reshape ( Hist_ab, [], 1 );

Idx3D = find ( ( Hist_3D ( :, 1 ) > 0 ) | ( Hist_3D ( :, 2 ) > 0 ) | ( Hist_3D ( :, 3 ) > 0 ) );
Hist_3D = Hist_3D ( Idx3D, : );

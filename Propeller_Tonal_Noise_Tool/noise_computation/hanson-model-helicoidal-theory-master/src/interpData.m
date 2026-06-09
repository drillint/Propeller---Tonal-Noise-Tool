%{
Script for interpolating input data while taking care of NaN and Inf values.

Copyright 2023 Jatinder Goyal

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
%}
%%
function outData = interpData(inData,xInterp)

    InfNaNidx = isinf(inData(:,2)) + isnan(inData(:,2));
    outData = interp1(inData(~InfNaNidx,1),inData(~InfNaNidx,2),xInterp,'spline',NaN); 
    NaNidx = isnan(outData);
    outData(NaNidx) = interp1(xInterp(~NaNidx),outData(~NaNidx),xInterp(NaNidx),'nearest','extrap');
%{
Main script for running the Hanson's tonal noise model in frequency domain.

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
clearvars
clc
UserInput
inputHanson
%% Harmonics
m = 1:6; %harmonic

%% calculate acoustic pressure
P = hanson(Amb,prop,op,obs,dep,m);

%{
%% post-processing
steps = 500;
t = timeSignal(prop,op,obs,m,P,steps); % calcuate the pressure time signal
figure;
plot(t.range,t.p) %plot pressure-time signal
% hold on
%}
SPL = SPLprocessing(Amb,P); %Calculate SPL and OASPL

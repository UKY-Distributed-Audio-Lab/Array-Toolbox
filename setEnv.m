function env = setEnv(c,mpos,fov);
% function setEnv returns data structure of environmental parameters
% the scripts that calls this function will have the env paramaters for the
% date of the recording
% input paramaters are setPro(c,mpos,fov)
% env - data structure describing the environmental parameters
%       fields include:
%           1) c - speed of sound
%           2) mpos - microphone positions
%           3) fov - field of view

env = struct('c',[c],'mpos',[mpos],'fov',[fov]);
%%Setup

% Clear the workspace
close all;
clearvars;
sca;

%Set to '2' for testing purposes. Change to '0' during trials (on windows
%PC)
Screen('Preference', 'SkipSyncTests', 2)

% Setup PTB with some default values
PsychDefaultSetup(2); 

%makes screen translucent
%PsychDebugWindowConfiguration( );

%getting ID of participant
ParticipantID = (cell2mat(inputdlg("Participant ID")));
while isempty(ParticipantID)
    ParticipantID = (cell2mat(inputdlg("You need to provide a participant ID:")));
end
%scales images
imgscale = 1;

%deactivate key input to editor
%ListenChar(-1);

%get folder where images are stored
imageFolder = '101 faces';

%retrieve image file names
imgArray = dir(fullfile(imageFolder,'*.jpg'));
imgList = {imgArray(:).name};

%getting number of trials = number of images in folder
nTrials = length(imgList);

%Stimulus values
imgArr = cell2mat(imgList');
preStimValues = str2num(imgArr(:,7:9))';
%standardizing stim values to 0-100.
stimValues = (preStimValues-min(preStimValues)) / (max(preStimValues)-min(preStimValues)) * 100;

% Randomize the trial list
randomizedTrials = randperm(nTrials);

% Set the screen number to the external secondary monitor if there is one
% connected
screenNumber = max(Screen('Screens'));
HideCursor(screenNumber);
% Define black, white and grey
white = WhiteIndex(screenNumber);
grey = white / 2;
black = BlackIndex(screenNumber);

% Open the screen
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, [] ,  32 ,  2,...
    [], [],  kPsychNeed32BPCFloat);
% Flip to clear
Screen('Flip', window);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);
%screen size in pixels
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Set the text size
Screen('TextSize', window, 40);

% Query the maximum priority level
topPriorityLevel = MaxPriority(window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);



% Now we set the number of times we want to do each condition, then make a
% full condition vector and then shuffle it. This will randomly order the
% orientation we present our Gabor with on each trial.
%numRepeats = 1;
%condVector = Shuffle(repmat(stimValues, 1, numRepeats));

% Calculate the number of trials
%numTrials = numel(condVector);

% Make a vector to record the response for each trial
respVector = zeros(1, nTrials);   

% Make a vector to count how many times we present each stimulus. This is a
% good check to make sure we have done things right and helps us when we
% input the data to plot anf fit our psychometric function
%countVector = zeros(1, nTrials);

%----------------------------------------------------------------------
%                       Timing Information
%----------------------------------------------------------------------

% Presentation Time for the Gabor in seconds and frames
presTimeSecs = 0.5;
presTimeFrames = round(presTimeSecs / ifi);

% Interstimulus interval time in seconds and frames
isiTimeSecs = 2;
isiTimeFrames = round(isiTimeSecs / ifi);

% Numer of frames to wait before re-drawing
waitframes = 1;


%----------------------------------------------------------------------
%                       Keyboard information
%----------------------------------------------------------------------

% Define the keyboard keys that are listened for. We will be using the left
% and right arrow keys as response keys for the task and the escape key as
% a exit/reset key

escapeKey = KbName('ESCAPE');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');


%----------------------------------------------------------------------
%                       Experimental loop
%----------------------------------------------------------------------
    % If this is the first trial we present a start screen and wait for a
    % key-press
DrawFormattedText(window, ...
    'In dieser Aufgabe erscheinen Gesichter mehrerer Personen. \n\nEntscheide, ob die Person wütend oder glücklich ist.\nIst sie glücklich, drücke die LINKE Pfeiltaste.\nIst sie wütend, drücke die RECHTE Pfeiltaste.\n\n -- Weiter mit Leertaste -- ',...
    'center', 'center', black);
Screen('Flip', window);
KbStrokeWait;
DrawFormattedText(window, ...
    'Antworte schnell und spontan.\n\nIn einigen Fällen wird nicht eindeutig sein, ob die Person wütend oder glücklich ist.\nWähle die Emotion, die am Ehesten zutrifft.\n\n -- Weiter mit Leertaste -- ',...
    'center', 'center', black);
Screen('Flip', window);
KbStrokeWait;
DrawFormattedText(window, ...
    '-- Drücke die Leertaste, um zu starten --',...
    'center', screenYpixels-100, black);
Screen('Flip', window);
KbStrokeWait;
% Animation loop: we loop for the total number of trials
Screen('TextSize', window, 50);
for trial = 1:nTrials

    %theAngle = (randomizedTrials(trial)-6) * 0.1;

    % Change the blend function to draw an antialiased fixation point
    % in the centre of the screen
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
 

    % Flip again to sync us to the vertical retrace at the same time as
    % drawing our fixation point
    Screen('DrawDots', window, [xCenter; yCenter], 10, black, [], 2);
    vbl = Screen('Flip', window);

    % Now we present the isi interval with fixation point minus one frame
    % because we presented the fixation point once already when getting a
    % time stamp
    for frame = 1:isiTimeFrames - 1

       
        %Draw the happy and angry indicators
        DrawFormattedText(window,'WÜTEND','right', screenYpixels-100, black, [], [], [], [], [], [100, 100, screenXpixels-100, screenYpixels-100]);
        DrawFormattedText(window,'GLÜCKLICH',100, screenYpixels-100, black);
        % Draw the fixation point
        Screen('DrawDots', window, [xCenter; yCenter], 10, black, [], 2);

        % Flip to the screen
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
    end
    respToBeMade = true;
    startIteration(trial) = tic;
    % Now we draw the Face
    for frame = 1:presTimeFrames
        
        file = imgList{randomizedTrials(trial)};
        fimg = imread(fullfile(imageFolder,file));
        imageTexture = Screen('MakeTexture', window, fimg);
        [s1, s2, s3] = size(fimg);  
        Screen('DrawTexture', window, imageTexture, [], [screenXpixels/2 - s2*imgscale/2 screenYpixels/2 - s1*imgscale/2  xCenter+s2*imgscale/2 yCenter+s1*imgscale/2], 0);
        % Set the right blend function for drawing the gabors
        Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');
        %Draw the happy and angry indicators
        DrawFormattedText(window,'WÜTEND','right', screenYpixels-100, black, [], [], [], [], [], [100, 100, screenXpixels-100, screenYpixels-100]);
        DrawFormattedText(window,'GLÜCKLICH',100, screenYpixels-100, black);
        % Flip to the screen
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);

        [keyIsDown,secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            ShowCursor;
            ListenChar(0);
            sca;
            return
        elseif keyCode(leftKey)
            respTime(randomizedTrials(trial)) = toc(startIteration(trial));
            response = 0;
            respToBeMade = false;
            break;
        elseif keyCode(rightKey)
            respTime(randomizedTrials(trial)) = toc(startIteration(trial));
            response = 1;
            respToBeMade = false;
            break;
        end

    end

  
    % Change the blend function to draw an antialiased fixation point
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

    % Draw the fixation point
    Screen('DrawDots', window, [xCenter; yCenter], 10, black, [], 2);
    %Draw the happy and angry indicators
    DrawFormattedText(window,'WÜTEND','right', screenYpixels-100, black, [], [], [], [], [], [100, 100, screenXpixels-100, screenYpixels-100]);
    DrawFormattedText(window,'GLÜCKLICH',100, screenYpixels-100, black);
    % Flip to the screen
    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);

    % Now we wait for a keyboard button signaling the observers response.
    % The left arrow key signals an "angry" response and the right arrow key
    % a "happy" response. You can also press escape if you want to exit the
    % programm
 
    waitingtime = tic;
    while respToBeMade
        if toc(waitingtime) >= 2
            Screen('DrawDots', window, [xCenter; yCenter], 10, black, [], 2);
            DrawFormattedText(window, '\n\n\n\n\n\nplease give a response', 'center', 'center', black);
            %Draw the happy and angry indicators
            DrawFormattedText(window,'WÜTEND','right', screenYpixels-100, black, [], [], [], [], [], [100, 100, screenXpixels-100, screenYpixels-100]);
            DrawFormattedText(window,'GLÜCKLICH',100, screenYpixels-100, black);
            Screen('Flip', window);
        end
        [keyIsDown,secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            ShowCursor;
            ListenChar(0);
            sca;
            return
        elseif keyCode(leftKey)
            respTime(randomizedTrials(trial)) = toc(startIteration(trial));
            response = 0;
            respToBeMade = false;
        elseif keyCode(rightKey)
            respTime(randomizedTrials(trial)) = toc(startIteration(trial));
            response = 1;
            respToBeMade = false;
        end
        
    end

    % Record the response
    respVector(trial) = response;

    % Add one to the counter for th bat stimulus
    %countVector(stimValues == theAngle) = countVector(stimValues == theAngle) + 1;
    disp("Number of trials remaining: ");
    nrem = nTrials - trial;
    disp (nrem);
end  
imgNames = string(imgArr(randomizedTrials, :));
data = table(imgNames, 100-stimValues(randomizedTrials)', respVector', respTime', 'VariableNames', {'image', 'angerIntensity', 'respondAngry', 'reactionTime'});

figure;
dataray = table2array(data);
%preparing data for plotting
values = str2double(dataray(:,2:3));

mResp = zeros(101,2);
for intensity = 1:101
 %   insy = (intensity - 6) * 0.1;
    index = abs(values(:,1) - (intensity-1)) < eps;
    %index = (values(:,1) == exp);
    mResp(102-intensity, :) = [intensity-1, mean(values(index, 2))];
end
    
plot(mResp(:,2), 'ro-', 'MarkerFaceColor', 'r');
axis([0 100 min(mResp(:, 2)) max(mResp(:, 2))]);
xlabel('Happy --- Angry');
ylabel('IsAngry');
title('Psychometric function');
writetable(data, ParticipantID)
% Clean up
ListenChar(0);
sca;
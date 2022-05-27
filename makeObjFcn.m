function ObjFcn = makeObjFcn(XTrain,YTrain,XValidation,YValidation)
%{
Define the objective function for optimization. This function performs the following steps:
1. Takes the values of the optimization variables as inputs. 
 bayesopt calls the objective function with the current values of the optimization variables in a table with each column name equal to the variable name. For example, the current value of the network section depth is optVars.SectionDepth.
2. Defines the network architecture and training options.
3. Trains and validates the network.
4. Saves the trained network, the validation error, and the training options to disk.
5. Returns the validation error and the file name of the saved network.
%}

imageSize = size(XTrain,1:3); %eg. [75,40,1] hxwxc
ObjFcn = @valErrorFun;
    function [valError,cons,fileName] = valErrorFun(optVars)
        [XTrain,YTrain] = AH_SMOTE(XTrain, YTrain, optVars.kSMOTE); % Use kNN to amplify dataset
        nClass = numel(unique(YTrain));
        depth = optVars.SectionDepth;
        %numF = round(16/sqrt(optVars.SectionDepth));
        if depth == 1
            layers = [
            imageInputLayer(imageSize)
            convolution2dLayer(round(optVars.filterSize),round(optVars.numFilters),'Padding',1) %first layer: 16 3x3 filters
            batchNormalizationLayer
            reluLayer
            maxPooling2dLayer(2,'Stride',2);
            
            fullyConnectedLayer(nClass)
            softmaxLayer
            classificationLayer];
        
        elseif depth == 2            
            layers = [
            imageInputLayer(imageSize)
            
            % The spatial input and output sizes of these convolutional
            % layers are 32-by-32, and the following max pooling layer
            % reduces this to 16-by-16.
            
            convolution2dLayer(round(optVars.filterSize),round(optVars.numFilters),'Padding',1) %first layer: 16 3x3 filters
            batchNormalizationLayer
            reluLayer
            maxPooling2dLayer(2,'Stride',2);
            %maxPooling2dLayer(3,'Stride',2,'Padding','same')
            
            % The spatial input and output sizes of these convolutional
            % layers are 16-by-16, and the following max pooling layer
            % reduces this to 8-by-8.
            convolution2dLayer(round(optVars.filterSize),round(2*optVars.numFilters),'Padding',1) %first layer: 16 3x3 filters
            batchNormalizationLayer
            reluLayer
            %maxPooling2dLayer(2,'Stride',2);
            %maxPooling2dLayer(3,'Stride',2,'Padding','same')
            
            % Add the fully connected layer and the final softmax and
            % classification layers.
            fullyConnectedLayer(nClass)
            softmaxLayer
            classificationLayer];
        
        elseif depth == 3            
            layers = [
            imageInputLayer(imageSize)           
            
            convolution2dLayer(round(optVars.filterSize),round(optVars.numFilters),'Padding',1) %first layer: 16 3x3 filters
            batchNormalizationLayer
            reluLayer
            maxPooling2dLayer(2,'Stride',2);
             
            convolution2dLayer(round(optVars.filterSize),round(2*optVars.numFilters),'Padding',1) 
            batchNormalizationLayer
            reluLayer
            maxPooling2dLayer(2,'Stride',2);
            
            convolution2dLayer(round(optVars.filterSize),round(2*optVars.numFilters),'Padding',1)
            batchNormalizationLayer
            reluLayer
            
            % Add the fully connected layer and the final softmax and
            % classification layers.
            fullyConnectedLayer(nClass)
            softmaxLayer
            classificationLayer];
        end
%{
Specify options for network training. Optimize the initial learning rate, SGD momentum, and L2 regularization strength.
Specify validation data and choose the 'ValidationFrequency' value such that trainNetwork validates the network once per epoch. Train for a fixed number of epochs and lower the learning rate by a factor of 10 during the last epochs. This reduces the noise of the parameter updates and lets the network parameters settle down closer to a minimum of the loss function.
%}
        miniBatchSize = 128;
        validationFrequency = floor(numel(YTrain)/miniBatchSize);
        options = trainingOptions('sgdm', ...
            'InitialLearnRate',optVars.InitialLearnRate, ... 
            'LearnRateSchedule','piecewise', ...
            'LearnRateDropPeriod',8, ...
            'LearnRateDropFactor',0.1, ...
            'MiniBatchSize',miniBatchSize, ...
            'L2Regularization',optVars.L2Regularization, ...  
            'MaxEpochs',30, ...
            'Verbose',false, ...
            'Plots','training-progress', ...
            'ValidationData',{XValidation,YValidation}, ...
            'ValidationFrequency',validationFrequency);
            %'Shuffle','every-epoch', ...
            % 'Momentum',optVars.Momentum, ...
       % Train the network and plot the training progress during training. Close all training plots after training finishes.
        trainedNet = trainNetwork(XTrain,YTrain,layers,options);
        close(findall(groot,'Tag','NNET_CNN_TRAININGPLOT_FIGURE'))
        % Evaluate the trained network on the validation set, calculate the predicted image labels, and calculate the error rate on the validation data.
        YPredicted = classify(trainedNet,XValidation);
        valError = 1 - mean(YPredicted == YValidation);

%{
Create a file name containing the validation error, and save the network, 
validation error, and training options to disk. 
The objective function returns fileName as an output argument, 
and bayesopt returns all the file names in BayesObject.UserDataTrace. 
The additional required output argument cons specifies constraints among the variables. There are no variable constraints.
%}
        fileName = "paraTuning_" + num2str(valError) + ".mat";
        save(fileName,'trainedNet','valError','options')
        cons = [];
        
    end
end


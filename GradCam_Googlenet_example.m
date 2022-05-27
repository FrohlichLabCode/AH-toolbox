% original example:
Gnet = googlenet;
inputSize = net.Layers(1).InputSize(1:2);
img = imread("sherlock.jpg");
img = imresize(img,inputSize);
[classfn,score] = classify(net,img);
imshow(img);
title(sprintf("%s (%.2f)", classfn, score(classfn)));
lgraph = layerGraph(Gnet);
lgraph = removeLayers(lgraph, lgraph.Layers(end).Name); % remove classification output layer
dlnet = dlnetwork(lgraph);
softmaxName = 'prob';
featureLayerName = 'inception_5b-output';
dlImg = dlarray(single(img),'SSC');
[featureMap, dScoresdMap] = dlfeval(@gradcam, dlnet, dlImg, softmaxName, featureLayerName, classfn);
gradcamMap = sum(featureMap .* sum(dScoresdMap, [1 2]), 3);
gradcamMap = extractdata(gradcamMap);
gradcamMap = rescale(gradcamMap);
gradcamMap = imresize(gradcamMap, inputSize, 'Method', 'bicubic');
imshow(img);
hold on;
imagesc(gradcamMap,'AlphaData',0.5);
colormap jet
hold off;
title("Grad-CAM");


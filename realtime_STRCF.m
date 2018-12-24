clear;
clc;
cam = webcam();
setup_paths();

count = 1;
time = 0;
while true
    tic();
    img = snapshot(cam);
    img = flip(img,2);
    img_gray = rgb2gray(img);
    %frameSize = size(img);
    % Create a cascade detector object.
    faceDetector = vision.CascadeObjectDetector();
    % Read a video frame and run the detector.
    bbox = step(faceDetector, img_gray);
    [seq] = create_realtime_sequence(img_gray,bbox,img);
    if bbox
        results = run_realtime_STRCF(seq,cam);
    else
        imagesc(img);
        hold on;
        text(10, 10, ['FPS: ' int2str(count/time)], 'color', [0 1 1]);
        hold off;
        axis off;axis image;set(gca, 'Units', 'normalized', 'Position', [0 0 1 1]);
    end
    time = time + toc();
    % Draw the returned bounding box around the detected face.
    %imgOut = insertObjectAnnotation(img,'rectangle',bbox,'Face');
    %imshow(imgOut);
end



% 
% % This demo script runs the STRCF tracker with deep features on the
% % included "Human3" video.
% 
% % Add paths
% setup_paths();
% 
% %  Load video information
% base_path  =  './sequences';
% %video  = choose_video(base_path);
% video = 'Human3';
% 
% video_path = [base_path '/' video];
% [seq, gt_boxes] = load_video_info(video_path);
% 
% % Run STRCF
% results = run_STRCF(seq);
% %results = run_DeepSTRCF(seq);
% 
% pd_boxes = results.res;
% thresholdSetOverlap = 0: 0.05 : 1;
% success_num_overlap = zeros(1, numel(thresholdSetOverlap));
% res = calcRectInt(gt_boxes, pd_boxes);
% 
% for t = 1: length(thresholdSetOverlap)
%     success_num_overlap(1, t) = sum(res > thresholdSetOverlap(t));
% end
% cur_AUC = mean(success_num_overlap) / size(gt_boxes, 1);
% FPS_vid = results.fps;
% display([video  '---->' '   FPS:   ' num2str(FPS_vid)   '    op:   '   num2str(cur_AUC)]);
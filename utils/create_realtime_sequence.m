function [params] = create_realtime_sequence(img,bbox,img_real)

    seq.format = 'realtime';
    seq.init_rect = bbox;
    seq.img = img;
    seq.num_frames = 1;
    seq.img_real = img_real;
    seq.isFirst = false;
    seq.isOldPos = false;
    params.seq = seq;

end


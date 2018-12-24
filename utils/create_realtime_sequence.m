function [params] = create_realtime_sequence(img,bbox)

    seq.format = 'realtime';
    seq.init_rect = bbox;
    seq.img = img;
    seq.num_frames = 1;
    seq.isFirst = false;
    seq.isOldPos = false;
    params.seq = seq;

end


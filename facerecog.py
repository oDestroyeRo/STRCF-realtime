import face_recognition
import cv2
import numpy as np
import PIL
import glob

def init():
    global known_face_encodings,known_face_names
    """Return list of words containing 'son'"""
    known_face_encodings = []
    known_face_names = []

    for filename in glob.glob('images/*.jpg'):
        image = face_recognition.load_image_file(filename)
        face_encoding = face_recognition.face_encodings(image)[0]
        f = filename.split(".")
        f = f[0].split("/")
        known_face_encodings.append(face_encoding)
        known_face_names.append(f[1])

def recog(frame):
    image = np.asarray(frame)
    small_frame = cv2.resize(image, (0, 0), fx=0.50, fy=0.50)
    rgb_small_frame = small_frame[:, :, ::-1]
    face_locations = face_recognition.face_locations(rgb_small_frame, model="cnn")
    face_encodings = face_recognition.face_encodings(rgb_small_frame, face_locations)
    face_names = []
    for face_encoding in face_encodings:
        matches = face_recognition.compare_faces(known_face_encodings, face_encoding)
        name = "Unknown"
        if True in matches:
            first_match_index = matches.index(True)
            name = known_face_names[first_match_index]

        face_names.append(name)
    return zip(face_locations, face_names)

def detect(frame):
    image = np.asarray(frame)
    small_frame = cv2.resize(image, (0, 0), fx=0.50, fy=0.50)
    rgb_small_frame = small_frame[:, :, ::-1]
    face_locations = face_recognition.face_locations(rgb_small_frame)
    return face_locations

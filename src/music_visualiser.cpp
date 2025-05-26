// text i/o
#include <iostream>
// audio i/o
#include "AudioFile.h"
// image i/o
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "EasyBMP.h"
//command line tools
#include <cstdlib>

// math
#include <cmath>
#include <complex>

// FFT
const double PI = acos(-1);

using namespace std;

template<typename T>
void fft(vector<complex<T>> & a, bool invert);

template<typename T>
void grab_audio_data(vector< vector<T> >& samples, int& num_channels, int& num_samples, int& sample_rate, const string filename);

string fixed_len(int x, int len) // only positive x, lol !
{
    const auto digits = int(log10(abs(x)));
    string str = "";
    for(int i=0; i<len-1-digits; i++) str +="0";
    str += to_string(x);
    //cout << x << " has " << digits << " digits so " << str << " has " << len << " digits." << endl;
    return str;
}
template<typename T>
void output_graph(vector<complex<T>> & values, int wdth, int hght, const string& frame_string);

void test_image();

void test_fft();

int main()
{
    test_image();

    /*
    string audio = "..\\wav_inputs\\a_lady.wav";
    string video = "..\\outputs\\a_lady.mp4";
    string video_with_sound = "..\\outputs\\a_lady_with_sound.mp4";
    */
    /*
    string audio = "..\\wav_inputs\\21_guns.wav";
    string video = "..\\outputs\\21_guns.mp4";
    string video_with_sound = "..\\outputs\\21_guns_with_sound.mp4";
    */
    ///*
    string audio = "..\\wav_inputs\\golden_brown.wav";
    string video = "..\\outputs\\golden_brown.mp4";
    string video_with_sound = "..\\outputs\\golden_brown_with_sound.mp4";
    //*/
    vector< vector<double> > samples; 
    int num_channels, num_samples, sample_rate;
    grab_audio_data(samples, num_channels, num_samples, sample_rate, audio);
    if (num_channels==1) samples.push_back(samples[0]); // If it's mono, just double the channel, lol
    double length_in_seconds = double(num_samples)/sample_rate;

    double frames_per_second = 12.0;
    double block_length_in_seconds = 1.0/frames_per_second;
    int block_length = block_length_in_seconds * sample_rate;
    cout << "Target FPS: " << frames_per_second << endl;
    cout << "Target block_length_in_seconds: " << block_length_in_seconds << endl;
    cout << "Target block_length: " << block_length << endl;
    block_length = pow(2, int(log2(block_length))); //rescale to an efficient power of 2
    block_length_in_seconds = block_length / double(sample_rate);
    frames_per_second = 1.0/block_length_in_seconds;
    cout << "Actual FPS: " << frames_per_second << endl;
    cout << "Actual block_length_in_seconds: " << block_length_in_seconds << endl;
    cout << "Actual block_length: " << block_length << endl;

    int num_blocks = int(num_samples / block_length); // How many images to generate
    cout << "Frame count: " << num_blocks << endl;

    vector<complex<double>> block_audio (block_length, 0); double block_max;

    int hght = 1000;
    int wdth = block_length/2;
    using Row = vector<uint8_t>;
    Row row(wdth * 3, 0); // Just for initialising
    vector<Row> ring_buffer (hght, row);
    int ring_head = 0; // points to the "first" row in the buffer

    vector<uint8_t> image_out(hght * wdth * 3);

    int str_len = log10(abs(num_blocks)) + 1;
    string out_str;
    for(int block=0; block<num_blocks; block++)
    {
        for(int i=0; i<block_length; i++)
            block_audio[i] = 1.0*samples[0][block*block_length+i] + 1.0i*samples[1][block*block_length+i]; //left channel real, right channel imaginary....

        fft(block_audio, false);

        block_max=0; 
        for (const auto & val : block_audio) if (abs(val)>block_max) block_max=abs(val); // find maximum
        if (block_max>0) for (auto & val : block_audio) val = 255.0*(abs(val)/block_max); // normalise it to 0.0-255.0
        for (int x=0; x<wdth; x++) // write line to the image buffer
        {
            ring_buffer[ring_head][3*x + 0] = block_audio[x].real(); // R
            ring_buffer[ring_head][3*x + 1] = block_audio[x].real(); // G
            ring_buffer[ring_head][3*x + 2] = block_audio[x].real(); // B
        }
        ring_head = (ring_head + 1)%hght; // advance the head

        // flatten buffer into output array
        for (int y=0; y<hght; y++) copy(ring_buffer[(ring_head+y)%hght].begin(), ring_buffer[(ring_head+y)%hght].end(), image_out.begin() + y * wdth * 3);

        // write frame PNG
        out_str = "..\\output_frames\\"+fixed_len(block, str_len)+"graph.png";
        stbi_write_png(out_str.c_str(), wdth, hght, 3, image_out.data(), wdth * 3);

        cout << "frame " << block << " done!" << endl;
    }

    cout << "Attempting ffmpeg..." << endl;

    string temp = "";
    int ffmpeg_fail_code;
    //Compile Frames with sequantial filenames <PRE><NUM><POST> (where NUM is a number of width N):
    temp = "ffmpeg -y -framerate 21.5 -i ..\\output_frames\\%0" + to_string(str_len) + "d" + "graph.png " + video + " 2>&1"; 
    ffmpeg_fail_code = system(temp.c_str()); 
    if (ffmpeg_fail_code) cerr << "ffmpeg video creation failed with code " << ffmpeg_fail_code << endl;
    //Add audio:
    temp = "ffmpeg -y -i " + video + " -i " + audio + " -map 0:v -map 1:a -c:v copy -shortest " + video_with_sound + " 2>&1";
    ffmpeg_fail_code = system(temp.c_str()); 
    if (ffmpeg_fail_code) cerr << "ffmpeg audio adding failed with code " << ffmpeg_fail_code << endl;

    cout << flush; cin.clear(); // reset in case of weird console state from system calls
    cout << "ffmpeg conversion done!" << endl;

    return 0;
}

template<typename T>
void grab_audio_data(vector< vector<T> >& samples, int& num_channels, int& num_samples, int& sample_rate, const string filename)
{
    // Create an AudioFile object
    AudioFile<T> audioFile;

    // Load an audio file
    bool loadedOK = audioFile.load (filename);
    assert (loadedOK);

    samples = audioFile.samples;
    //audioFile.samples is vector<vector<type>> of shape (num_channels, num_samples)
    //audioFile.samples[channel][sampleIndex] is <type>

    // Get some information about the loaded audio
    num_channels = audioFile.getNumChannels(); // number of channels
    num_samples = audioFile.getNumSamplesPerChannel(); // total number of samples
    sample_rate = audioFile.getSampleRate(); // samples per second
    cout << "audioFile.getLengthInSeconds() = " << audioFile.getLengthInSeconds() << " =?= " << double(num_samples)/sample_rate << " = num_samples/sample_rate" << endl;
    cout << "audioFile.getBitDepth() = " << audioFile.getBitDepth() << " <?= " << 32 << " = float bit depth" << endl;

    // or, just use this quick shortcut to print a summary to the console
    audioFile.printSummary();

    return;
}

template<typename T>
void fft(vector<complex<T>> & a, bool invert) {
    int n = a.size();

    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;

        if (i < j)
            swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? -1 : 1);
        complex<T> wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            complex<T> w(1);
            for (int j = 0; j < len / 2; j++) {
                complex<T> u = a[i+j], v = a[i+j+len/2] * w;
                a[i+j] = u + v;
                a[i+j+len/2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert) {
        for (complex<T> & x : a)
            x /= n;
    }
}

template<typename T>
void output_graph(vector<complex<T>> & values, int wdth, int hght, const string& frame_str)
{
    // Parameters
    int bit_depth = 32;
    // Declare a new bitmap object
    BMP image;
    // Set size to 640 � 480
    image.SetSize(wdth, hght);
    // Set its color depth to 8-bits
    image.SetBitDepth(bit_depth);

    // Change colors
    for(int x=0; x<wdth; x++)
        for(int y=0; y<hght; y++)
            image.SetPixel(x, y, {0, 0, 0, 255});

    int N = values.size();

    double max_val = 0;
    for(int i=0; i<N; i+=1)
        if(abs(values[i])>max_val)
            max_val = abs(values[i]);
    //if(max_val == 0.0) max_val = 1.0;

    int i, yr, yi;
    for(int x=0; x<wdth; x++)
    {
        if (max_val > 0) // most cases
        {
            i = 0.5*double(N)*x/wdth; //transform from image scale to vector index (halved because symmetric)
            yr = int(0.5*hght + 0.5*hght*values[i].real()/max_val);
            yi = int(0.5*hght + 0.5*hght*values[i].imag()/max_val);
            if (abs(yr)>=hght) yr=int(((yr>0)-(yr<0)+1)*hght/2.0); // If I've fucked up and overflowed somewhere
            if (abs(yi)>=hght) yi=int(((yi>0)-(yi<0)+1)*hght/2.0); // If I've fucked up and overflowed somewhere
        }
        else // sometimes there's literally no sound
        {
            yr = 0.0;
            yi = 0.0;
        }

        image.SetPixel(x, yr, {0, 0, 255, 255}); //red for real
        image.SetPixel(x, yi, {255, 0, 0, 255}); //blue for imaginary

        for (int h=0; h<abs(yr-hght/2); h++)
            image.SetPixel(x, int(0.5*hght+((yr-hght/2>0)-(yr-hght/2<0))*h), {0, 0, 255, 255}); //draw red line
        for (int h=0; h<abs(yi-hght/2); h++) 
            image.SetPixel(x, int(0.5*hght+((yi-hght/2>0)-(yi-hght/2<0))*h), {255, 0, 0, 255}); //draw blue line
    }

    // Save to file
    //cout << "wrote " << frame_str << endl;
    image.WriteToFile(frame_str.c_str());
}

void set_buffer_pixel(vector<uint8_t> & image_buffer, int wdth, int x, int y, uint8_t r, uint8_t g, uint8_t b)
{
    //Manual indexing :(
    int idx = (y*wdth + x) * 3;
    image_buffer[idx + 0] = r;
    image_buffer[idx + 1] = g;
    image_buffer[idx + 2] = b;
    return;
}
void test_image()
{
    cout << "generating test image" << endl;
    // Parameters
    int wdth = 1024;
    int hght = 1024;
    // Declare a new buffer for a full image
    vector<uint8_t> image_buffer(wdth * hght * 3); // x * y * rgb // alternately do wdth*hght*3 for rgba
    // Change colors
    uint8_t r, g, b;
    for(int x=0; x<wdth; x++)
    {
        for(int y=0; y<hght; y++)
        {
            r = (255.0*x)/wdth; g = 0; b = (255.0*y)/hght;
            //cout << r << " ";
            set_buffer_pixel(image_buffer, wdth, x, y, r, g, b);
        } //cout << endl;
    }
    // Save to file
    stbi_write_png("..\\outputs\\test.png", wdth, hght, 3, image_buffer.data(), wdth * 3);

    cout << "Image test done." << endl;

    return;
}

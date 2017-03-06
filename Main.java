import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.io.File;
import java.util.Arrays;

public class Main {
    public static BufferedImage loadImage(String fileName) {
        BufferedImage img = null;
        try {
            img = ImageIO.read(new File(fileName));
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
        return img;
    }
    
    public static void writeImage(BufferedImage img, String fileName) {
         try {
            ImageIO.write(img, "png", new File(fileName));
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
    
    public static BufferedImage processPixelwise(BufferedImage in, PixelTransformFunction f) {
        BufferedImage out = new BufferedImage(in.getWidth(),
                in.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
        WritableRaster inRaster = in.getRaster();
        WritableRaster outRaster = out.getRaster();
        for (int i = 0; i < in.getWidth(); i++) {
            for (int j = 0; j < in.getHeight(); j++) {
                int v = inRaster.getSample(i, j, 0);
                v = (int) (255 * f.transform(v / 255.0));
                outRaster.setSample(i, j, 0, v);
            }
        }
        return out;
    }
    
    public static BufferedImage processPixelwiseExponent(BufferedImage in, double exponent) {
        PixelTransformFunction transform = new PixelTransformFunction() {
            @Override
            public double transform(double v) {
                return Math.pow(v, exponent);
            }
        };
        return processPixelwise(in, transform);
    }
    
    public static double computeMaskAtPixel(
            WritableRaster raster, int x, int y, double[][] mask) {
        int maskOffset = (mask.length - 1) / 2;
        int xStart = x - maskOffset, yStart = y - maskOffset;
        if (xStart < 0 || yStart < 0 
                || xStart + mask.length > raster.getWidth()
                || yStart + mask.length > raster.getHeight())
            return 0.0;
        double sum = 0.0;
        for (int i = 0; i < mask.length; i++) {
            for (int j = 0; j < mask.length; j++) {
                sum += raster.getSample(xStart + i, yStart + j, 0) * mask[i][j];
            }
        }
        return sum;
    }
    
    public static BufferedImage processMask(BufferedImage in, double[][] mask) {
        BufferedImage out = new BufferedImage(in.getWidth(),
                in.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
        WritableRaster inRaster = in.getRaster();
        WritableRaster outRaster = out.getRaster();
        for (int i = 0; i < in.getWidth(); i++) {
            for (int j = 0; j < in.getHeight(); j++) {
                int v = (int) computeMaskAtPixel(inRaster, i, j, mask);
                if (v < 0) v = 0;
                if (v > 255) v = 255;
                outRaster.setSample(i, j, 0, v);
            }
        }
        return out;
    }
    
    public static double computeMedianAtPixel(
            WritableRaster raster, int x, int y, int size) {
        int offset = (size - 1) / 2;
        int xStart = x - offset, yStart = y - offset;
        if (xStart < 0 || yStart < 0 
                || xStart + size > raster.getWidth()
                || yStart + size > raster.getHeight())
            return 0.0;
        double[] values = new double[size * size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                values[i * size + j] = raster.getSample(xStart + i, yStart + j, 0);
            }
        }
        Arrays.sort(values);
        return values[size * size / 2 - 1];
    }
    
    public static BufferedImage processMedian(BufferedImage in, int size) {
        BufferedImage out = new BufferedImage(in.getWidth(),
                in.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
        WritableRaster inRaster = in.getRaster();
        WritableRaster outRaster = out.getRaster();
        for (int i = 0; i < in.getWidth(); i++) {
            for (int j = 0; j < in.getHeight(); j++) {
                int v = (int) computeMedianAtPixel(inRaster, i, j, size);
                outRaster.setSample(i, j, 0, v);
            }
        }
        return out;
    }
    
    public static double[][] lineMask(int size, int nLines) {
        double v = 1.0 / (size * nLines);
        double[][] mask = new double[size][size];
        for (int i = 0; i < nLines; i++) {
            for (int j = 0; j < size; j++) {
                mask[i][j] = v;
            }
        }
        return mask;
    }
    
    public static double[] computeDFTAtPixel(WritableRaster raster, int x, int y) {
        double re = 0.0, im = 0.0;
        for (int i = 0; i < raster.getWidth(); i++) {
            for (int j = 0; j < raster.getHeight(); j++) {
                double k = i * x / ((double) raster.getWidth());
                k += j * y / ((double) raster.getHeight());
                //re += Math.cos(-2 * Math.PI * k) * raster.getSample(i, j, 0);
                //im += Math.sin(-2 * Math.PI * k) * raster.getSample(i, j, 0);
                re += Math.cos(-2 * Math.PI * k) * raster.getSample(i, j, 0) * Math.pow(-1.0, i + j);
                im += Math.sin(-2 * Math.PI * k) * raster.getSample(i, j, 0) * Math.pow(-1.0, i + j);
            }
        }
        re /= raster.getWidth() * raster.getHeight();
        im /= raster.getWidth() * raster.getHeight();
        return new double[] {re, im};
    }
    
    public static double[][] computeDFT(BufferedImage img) {
        double[] re = new double[img.getWidth() * img.getHeight()];
        double[] im = new double[img.getWidth() * img.getHeight()];
        for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getHeight(); j++) {
                double[] dftValue = computeDFTAtPixel(img.getRaster(), i, j);
                re[i * img.getWidth() + j] = dftValue[0];
                im[i * img.getWidth() + j] = dftValue[1];
            }
            System.out.println(i);
        }
        return new double[][] {re, im};
    }
    
    public static BufferedImage makeImage(double[] data, double min, double max) {
        int size = (int) Math.sqrt(data.length);
        BufferedImage img = new BufferedImage(size, size,
                BufferedImage.TYPE_BYTE_GRAY);
        WritableRaster raster = img.getRaster();
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                int v = (int) (255 * ((data[i * size + j] - min) / (max - min)));
                if (v < 0) v = 0;
                if (v > 255) v = 255;
                raster.setSample(i, j, 0, v);
            }
        }
        return img;
    }
    
    public static double[] absValues(double[][] data) {
        double[] out = new double[data[0].length];
        for (int i = 0; i < data[0].length; i++) {
            out[i] = Math.sqrt(data[0][i] * data[0][i] + data[1][i] * data[1][i]);
        }
        return out;
    }
    
    public static void main(String[] args) {
        BufferedImage img = loadImage(args[0]);
        //img = processPixelwiseExponent(img, Double.parseDouble(args[2]));
        //img = processMask(img, new double[][] {{1, 1, 1}, {1, -8, 1}, {1, 1, 1}});
        //img = processMask(img, constantMask(Integer.parseInt(args[2])));
        //img = processMask(img, lineMask(21, 5));
        //img = processMedian(img, 3);
        //img = processMask(img, lineMask(3, 3));
        //writeImage(img, args[1]);
        double[][] dft = computeDFT(img);
        double[] abs = absValues(dft);
        writeImage(makeImage(abs, 0.0, 1.0), args[1]);
    }
}

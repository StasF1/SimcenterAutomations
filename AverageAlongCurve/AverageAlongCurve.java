/*
    API:            Simcenter STAR-CCM+ 15.04.010
    Project:        https://github.com/StasF1/StarCcmMacros
    License:        GNU General Public License 3.0 ( see LICENSE )
    Author:         Stanislav Stashevskii

    Macro:          AverageAlongCurve.java
    Description:    Avarage parameters along a curve
*/
package macro;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import star.base.neo.*;
import star.base.report.*;
import star.common.*;
import star.meshing.*;
import star.vis.*;

public class AverageAlongCurve extends StarMacro {
  public void execute() {
    Simulation simulation = getActiveSimulation();

    String regionOfPipe = "Assembly 1.big_coll";

    String[] fieldNames = new String[] { //FIXME : Post-process only requried field names
      "AbsoluteTotalPressure",
      "AbsolutePressure",
      "Density",
      "DynamicViscosity",
      "MassFlux",
      "Temperature",
      "TotalTemperature"
    };

    String pathToSaveCsv = MkdirFromSimulationName(simulation.getSessionPath(), ".PipeCuts");

    {
      List<double[]> origins = new ArrayList<double[]>();
        origins.add(new double[] {0.9, 2.9, 6.00});
        origins.add(new double[] {0.9, 2.9, 7.00});
        origins.add(new double[] {0.9, 2.9, 8.00});
        origins.add(new double[] {0.9, 2.9, 9.00});
        origins.add(new double[] {0.9, 2.9, 10.00});
        origins.add(new double[] {0.9, 2.9, 11.00});
        origins.add(new double[] {0.9, 2.9, 12.00});
        origins.add(new double[] {0.9, 2.9, 13.00});
        origins.add(new double[] {0.9, 2.9, 14.00});
        origins.add(new double[] {0.9, 2.9, 15.00});

      List<double[]> orientations = new ArrayList<double[]>();
        orientations.add(new double[] {0.0, 0, 1.0});
        orientations.add(new double[] {0.0, 0, 1.0});
        orientations.add(new double[] {0.0, 0, 1.0});
        orientations.add(new double[] {0.0, 0, 1.0});
        orientations.add(new double[] {0.0, 0, 1.0});
        orientations.add(new double[] {0.0, 0, 1.0});
        orientations.add(new double[] {0.0, 0, 1.0});
        orientations.add(new double[] {0.0, 0, 1.0});
        orientations.add(new double[] {0.0, 0, 1.0});
        orientations.add(new double[] {0.0, 0, 1.0});

      CreatePlaneSection(simulation, "alongCurveCut", regionOfPipe);
      CreateSurfaceAverageReport(simulation, "surfaceAverageAlongCurveCut", "alongCurveCut");

      for (String field : fieldNames) {
        String pipeCutsCsv = ConvertPipeCutsToCsv(
          CreatePipeCuts(simulation, regionOfPipe, field, origins, orientations)
        );
        SaveTextToFile(pathToSaveCsv + "\\" + field + ".csv", pipeCutsCsv);

        simulation.println(field + " CSV field by the tube length:\n" + pipeCutsCsv);
      }
      simulation.println("End");
    }
  }


  static String MkdirFromSimulationName(String sessionPath) {
    return MkdirFromSimulationName(sessionPath, "");
  }
  static String MkdirFromSimulationName(String sessionPath, String extension) {
    String sessionPathWoExtension = sessionPath.substring(0, sessionPath.lastIndexOf('.'));

    File pathToCreate = new File(sessionPathWoExtension + extension);
    boolean dirCreated = pathToCreate.mkdir();

    return sessionPathWoExtension + extension;
  }

  static void CreatePlaneSection(Simulation simulation,
                                 String presentationName, String regionName) {
    if (!simulation.getPartManager().has(presentationName)) {
      PlaneSection planeSection = (PlaneSection) simulation.getPartManager().createImplicitPart(
        new NeoObjectVector(new Object[] {}),
        new DoubleVector(new double[] {0.0, 0.0, 1.0}),
        new DoubleVector(new double[] {0.0, 0.0, 0.0}),
        0,
        1,
        new DoubleVector(new double[] {0.0})
      );

      Units units = ((Units) simulation.getUnitsManager().getObject("m"));

      planeSection.getOriginCoordinate().setCoordinate(
        units, units, units,
        new DoubleVector(new double[] {0.0, 0.0, 0.0})
      );

      planeSection.getInputParts().setQuery(null);
      Region region = simulation.getRegionManager().getRegion(regionName);
      planeSection.getInputParts().setObjects(region);
      planeSection.setPresentationName(presentationName);

      simulation.println("Created " + presentationName
                         + " plane section in the region " + regionName);
    }
  }

  static void CreateSurfaceAverageReport(Simulation simulation,
                                         String presentationName, String planeSectionName) {
    if (!simulation.getReportManager().has(presentationName)) {
      PlaneSection planeSection = 
        ((PlaneSection) simulation.getPartManager().getObject(planeSectionName));

      AreaAverageReport areaAverageReport = 
        simulation.getReportManager().createReport(AreaAverageReport.class);
      areaAverageReport.getParts().setQuery(null);
      areaAverageReport.getParts().setObjects(planeSection);
      areaAverageReport.setPresentationName(presentationName);

      simulation.println("Created " + presentationName
                         + " surface average report on the plane section " + planeSectionName);
    }
  }

  static void SaveTextToFile(String filename, String text) {
    // simulation.getSessionDirFile();
    try (PrintWriter out = new PrintWriter(filename)) {
      out.println(text);
    } catch (IOException exception) {
      System.out.println("Path do not exist --> [" + filename + "]");
    }
  }


  private class PipeCut {
    public PipeCut(double[] origin, double value) {
      coordinates_ = origin;
      value_ = value;
    }

    public double[] getCoordinates() {
      return coordinates_;
    }
    public double getValue() {
      return value_;
    }

    private double[] coordinates_;
    private double value_;
  };

  private void EditPlaneSection(Simulation simulation,
                                String regionOfPipe, double[] origin, double[] orientation) {
    PlaneSection planeSection =
      ((PlaneSection) simulation.getPartManager().getObject("alongCurveCut"));

    Units units = ((Units) simulation.getUnitsManager().getObject("m"));

    planeSection.getOriginCoordinate().setCoordinate(
      units, units, units,
      new DoubleVector(origin)
    );
    planeSection.getOrientationCoordinate().setCoordinate(
      units, units, units,
      new DoubleVector(orientation)
    );

    planeSection.getInputParts().setQuery(null);
    planeSection.getInputParts().setObjects(
      simulation.getRegionManager().getRegion(regionOfPipe)
    );
  }

  private double GetReportValue(Simulation simulation, String reportName, String fieldName) {
    AreaAverageReport surfaceAveragePipeCutReport =
      ((AreaAverageReport) simulation.getReportManager().getReport(reportName));

    PrimitiveFieldFunction field =
      ((PrimitiveFieldFunction) simulation.getFieldFunctionManager().getFunction(fieldName));

    surfaceAveragePipeCutReport.setFieldFunction(field);
    return surfaceAveragePipeCutReport.getValue();
  }

  private List<PipeCut> CreatePipeCuts(Simulation simulation,
                                       String regionOfPipe, String fieldName,
                                       List<double[]> origins, List<double[]> orientations) {
    List<PipeCut> pipeCuts = new ArrayList<PipeCut>();

    for (int i = 0; i < orientations.size(); i++) {
      double[] origin = origins.get(i);
      double[] orientation = orientations.get(i);

      EditPlaneSection(simulation, regionOfPipe, origin, orientation);
      pipeCuts.add(new PipeCut(
        origin,
        GetReportValue(simulation, "surfaceAverageAlongCurveCut", fieldName)
      ));
    }
    return pipeCuts;
  }

  private String ConvertPipeCutsToCsv(List<PipeCut> pipeCuts) {
    String pipeCutsCsv = new String("x,y,z,value\n");

    for (int i = 0; i < pipeCuts.size(); i++) {
      PipeCut pipeCut = pipeCuts.get(i);

      for (double coordinate : pipeCut.getCoordinates()) {
        pipeCutsCsv += String.valueOf(coordinate) + ",";
      }
      pipeCutsCsv += pipeCut.getValue() + "\n";
    }
    return pipeCutsCsv;
  }
}

/* (C) 2021 Stanislav Stasheuskii */
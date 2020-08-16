/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package multigeneediting;

import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JTable;

/**
 *
 * @author Yi
 */
public class Tool {
    public static String projectDir = new File("").getAbsolutePath();
    /*
    * 以行为单位读取文件，常用于读面向行的格式化文件
    */
    public static ArrayList<String> readFileByLines(String fileName) {
        File file = new File(fileName);
        BufferedReader reader = null;
        ArrayList<String> result=new ArrayList<String>();
        try {
            reader = new BufferedReader(new FileReader(file));
            String tempString = null;
            while ((tempString = reader.readLine()) != null) {
                result.add(tempString);
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e1) {
                }
            }
        }
        return result;
    }
    
    /*
    * Read file as array list
    */
    public static ArrayList<String[]> readFileAsArrayList(File file) {
        BufferedReader reader = null;
        ArrayList<String[]> result=new ArrayList();
        try {
            reader = new BufferedReader(new FileReader(file));
            String tempString = null;
            while((tempString = reader.readLine()) != null) {
                tempString = tempString.replaceAll(" ", "");
                if(tempString.length()>0){
                    String[]items = tempString.split(",");
                    result.add(items);
                }
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e1) {
                }
            }
        }
        result.remove(0);
        return result;
    }
    
    /*Call the python via command line*/
    public static void callPyByCommand(String commandStr){
        Process p = null; 
        try {
            p = Runtime.getRuntime().exec(commandStr);
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));   
            String readLine = null;   
            readLine = br.readLine();
            while (readLine != null) { 
                readLine = br.readLine();
                System.out.println(readLine); 
            } 
            if(br!=null){ 
                br.close();
            }
        } catch (IOException ex) {
            Logger.getLogger(Tool.class.getName()).log(Level.SEVERE, null, ex);
        }
        p.destroy(); 
        p=null; 
    }
    
    /*Copy the content*/
    public static void copyTableCell(JTable table){
        int col = table.getSelectedColumn();
        int row = table.getSelectedRow();
        if(col>-1&&row>-1){
            String selectedContent = table.getValueAt(row, col).toString().replaceAll("\n\n", "\n");
            Clipboard clip = Toolkit.getDefaultToolkit().getSystemClipboard();  
            Transferable tText = new StringSelection(selectedContent);  
            clip.setContents(tText, null);
        }
    }
}

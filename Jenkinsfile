import groovy.json.JsonSlurper
import hudson.util.RemotingDiagnostics
import jenkins.model.Jenkins
import java.util.regex.Pattern
import java.util.regex.Matcher

@groovy.transform.Field Properties config

// https://qa.nuxeo.org/jenkins/pipeline-syntax/globals

def PROJECT_NAME = "DPLASMA"

//node {
//    sh 'env > env.txt'
//    for (String i : readFile('env.txt').split("\r?\n")) {
//        println i
//    }
//}

config = new Properties()
//String configFileName = System.getenv('JENKINS_HOME') + "/dplasma.jenkins.conf"
// For some obscure reason the JENKINS_HOME is always NULL
String configFileName = "${env.JENKINS_HOME}/dplasma.jenkins.conf"
try {
    File configFile = new File(configFileName)
    config.load(configFile.newDataInputStream())
} catch (Exception ex) {
    println "Exception: " + ex.getMessage()
    error("Cant load local configuration file ${configFileName} (required for connection to bitbucket)")
}

propertiesData = [disableConcurrentBuilds()]
properties(propertiesData)

node {
    def regex = "^https://bitbucket\\.org/(.+?)/(.+?)/pull-requests/(\\d+)"
    try {
        def l = Pattern.compile(regex).matcher(env.CHANGE_URL)
        if( l.find() ) {
            config.put('organization', l.group(1))
            config.put('repository', l.group(2))
        }
        println "${config.organization} -> ${config.repository}"
    } catch (Exception ex) {
        println "Exception: " + ex.getMessage() + " (${env.CHANGE_URL})"
        println "Cant identify the organization and repository from ${env.CHANGE_URL} using ${regex}"
    }
    if( (null == config.get('organization')) || (null == config.get('repository')) ) {
        println "Missing organization or repository"
        error("Missing organization or repository")
    }
    if( null == config.get('useSlack') ) {
        config.put('useSlack', false)
    }
    String userpass = config.userName + ":" + config.userPassword;
    String basicAuth = "Basic " + userpass.bytes.encodeBase64().toString()
    config.put('basicAuth', basicAuth)
}

def COLOR_MAP = [
    'SUCCESS': 'good',
    'FAILURE': 'danger',
]
def slackResponse = null

// Save the master of the lack check
//git log --format="%H" -n 1 master

//currentBuild.properties.each { println "currentBuild.${it.key} -> ${it.value}" }
//propertiesData.properties.each { println "$it.key -> $it.value" }
//config.each { println "config.${it.key} -> ${it.value}" }

// To mark the starting of a new build the previous approval should be removed
node {
    unapprovePullRequest(config.repository, env.CHANGE_ID)
}

pipeline {
    agent any
    stages {
        stage ('Clone') {
            when {
                beforeAgent true  // execute the when clause early
                anyOf {
                    expression {
                        // https://github.com/jenkinsci/jenkins/blob/master/core/src/main/java/hudson/model/Result.java
                        return isApprovedBranch(config.repository, env.CHANGE_ID, env.CHANGE_AUTHOR)
                    }
                }
            }
            steps {
                script {
                    try {
                        blocks = [
                          [
                            "type": "section",
                            "text": [
                                "type": "mrkdwn",
                                "text": "${PROJECT_NAME} Jenkins: [${env.BUILD_URL}](${env.JOB_NAME} [${env.BUILD_NUMBER})"
                                    ]
                          ]]
                        slackResponse = slackSend(channel: '#ci', blocks:blocks)
                    } catch (Exception ex) {  //disable Slack
                        config.useSlack = false
                    }
                }
                checkout scm
            }
        }
        stage ('Fetch Submodule PaRSEC') {
            steps {
                sh "git --work-tree=${env.WORKSPACE} --git-dir=${env.WORKSPACE}/.git submodule update --force --init --recursive --remote"
            }
        }
        stage ('Build and Test') {
            environment {
                BUILDDIR = "build"
            }
            parallel {
               stage ('Debug w/mod-parsec') {
                    environment {
                        BUILDTYPE = "Debug"
                        BUILDDIR = "${env.BUILDDIR}/${env.BUILDTYPE}"
                    }
                    options {
                        timeout(time:1, unit: 'HOURS')
                    }
                    steps {
                        sh "echo 'Prepare for a ${env.BUILDTYPE} build'"
                        sh "mkdir -p ${env.BUILDDIR}"
                        dir (env.BUILDDIR) {
                            sh "${env.WORKSPACE}/contrib/jenkins/script.saturn ${env.BUILDTYPE}"
                        }
                    }
                }
                stage ('Release w/mod-parsec') {
                    environment {
                        BUILDTYPE = "Release"
                        BUILDDIR = "${env.BUILDDIR}/${env.BUILDTYPE}"
                    }
                    options {
                        timeout(time:1, unit: 'HOURS')
                    }
                    steps {
                        sh "echo 'Prepare for a ${env.BUILDTYPE} build'"
                        sh "mkdir -p ${env.BUILDDIR}"
                        dir (env.BUILDDIR) {
                            sh "${env.WORKSPACE}/contrib/jenkins/script.saturn ${env.BUILDTYPE}"
                        }
                    }
                }
                stage ('Debug w/ext-parsec') {
                    environment {
                        BUILDTYPE = "Debug-ext"
                        BUILDDIR = "${env.BUILDDIR}/${env.BUILDTYPE}"
                    }
                    options {
                        timeout(time:1, unit: 'HOURS')
                    }
                    steps {
                        sh "echo 'Prepare for a ${env.BUILDTYPE} w/external PaRSEC build'"
                        sh "mkdir -p ${env.BUILDDIR}"
                        dir (env.BUILDDIR) {
                            sh "${env.WORKSPACE}/contrib/jenkins/script.saturn ${env.BUILDTYPE}"
                        }
                    }
                }
            }
        }
    }
    post {
        always {
            script {
                BUILD_USER = env.CHANGE_AUTHOR
            }
            slackSend channel: '#ci',
                      color: COLOR_MAP[currentBuild.currentResult],
                      message: "*${currentBuild.currentResult}:* Job ${env.JOB_NAME} build ${env.BUILD_NUMBER} by ${BUILD_USER}\n More info at: ${env.BUILD_URL}"
        }
        regression {
            slackSend channel: slackResponse.channelId, timestamp: slackResponse.ts,
                      message: "${PROJECT_NAME} REGRESSION: Job <a href=\"${env.BUILD_URL}\">${env.JOB_NAME} [${env.BUILD_NUMBER}]</a><br>" +
                               "Pull Request <a href=\"${env.CHANGE_URL}\">#${env.CHANGE_ID}</a>."
        }
        success {
            approvePullRequest(config.repository, env.CHANGE_ID)
        }
        failure {
            unapprovePullRequest(config.repository, env.CHANGE_ID)
        }
    }
}

// Find the DPLASMA version from the CMakeLists.txt file
String getVersion() {
    try {
        def p = ['grep', '(DPLASMA_VERSION_', 'CMakeLists.txt'].execute() | ['awk', '/^set/ {printf \"%d.\", $3}'].execute()
        p.waitFor()
        if( 0 == p.exitValue() ) {
            return p.in.text[0..-2]
        } else {
            println "ERROR: ${p.exitValue()} -> ${p.err.text}"
        }
    } catch (err) {
        println "ERROR: while retrieving DPLASMA version from CMakeLists.txt"
        throw err
    }
    return "0.0"
}

def displayServerResponse( InputStream is ) {
    BufferedReader input = new BufferedReader(new InputStreamReader(is))
    String inputLine;
    while ((inputLine = input.readLine()) != null)
        println inputLine
    input.close();
}

// Check if the branch author has validated his work by approving
// the branch.
def isApprovedBranch(repository, pr, byWho) {
    Boolean prApproved = false

    url = "https://api.bitbucket.org/2.0/repositories/${config.organization}/${repository}/pullrequests/${pr}/"
    println "Check URL: ${url}"
    def url = new URL(url)

    def conn = url.openConnection()
    conn.setRequestProperty( "Authorization", config.basicAuth)

    try {
        statusCode = conn.getResponseCode()
        if (statusCode==200) {
            InputStream inputStream = conn.getInputStream()
            def names = new groovy.json.JsonSlurper().parseText(inputStream.text)
            //names.each{ k, v ->
            //    println "element name: ${k}, element value: ${v}"
            //}
            inputStream.close()
            names['participants'].each {
                if( it['approved'] ) {
                    if( it['user']['username'] == byWho ) {
                        prApproved = true
                        println "Pullrequest ${pr} approved by ${byWho}"
                    }
                }
            }
        } else {
            println "[Check Approval] Connection status code: $statusCode "
            println "[Check Approval] URL ${url.toString()}"
            println "[Check Approval] Server response:"
            println "[Check Approval] -----"
            response=displayServerResponse(conn.getErrorStream())
            println "[Check Approval] -----"
        }
    } catch (Exception e) {
        println "[Check Approval] Error connecting to the URL"
        println e.getMessage()
    } finally {
        if (conn != null) {
            conn.disconnect();
        }
    }
    return prApproved
}

def approvePullRequest(repository, pr) {
    Boolean prApproved = false

    def url = new URL("https://api.bitbucket.org/2.0/repositories/${config.organization}/${repository}/pullrequests/${pr}/approve")

    def HttpURLConnection conn = (HttpURLConnection)url.openConnection()
    conn.setRequestMethod("POST")
    conn.setRequestProperty( "Authorization", config.basicAuth)
    conn.setRequestProperty( "Content-type", "application/x-www-form-urlencoded")
    conn.setRequestProperty( "Content-Length", "0")
    conn.setDoOutput(true)
    conn.connect()

    try {
        statusCode = conn.getResponseCode()
        if (statusCode==200) {
            InputStream inputStream = conn.getInputStream()
            def names = new groovy.json.JsonSlurper().parseText(inputStream.text)
            //names.each{ k, v ->
            //    println "element name: ${k}, element value: ${v}"
            //}
            inputStream.close()
            prApproved = true
            println "[Approve RP] PR succesfully approved"
        } else if (statusCode==409) {
            println "[Approve RP] PR already approved by this user"
            prApproved = true
        } else {
            println "[Approve PR] Connection status code: $statusCode "
            println "[Approve PR] URL ${url.toString()}"
            println "[Approve PR] Server response:"
            response = displayServerResponse(conn.getErrorStream())
            println "[Approve PR] ${response}"
            println "[Approve PR] -----"
        }
    } catch (Exception e) {
        println "[Approve PR] Error connecting to the URL"
        println e.getMessage()
    } finally {
        if (conn != null) {
            conn.disconnect();
        }
    }
    return prApproved
}

def unapprovePullRequest(repository, pr) {
    Boolean prUnapproved = false

    def url = new URL("https://api.bitbucket.org/2.0/repositories/${config.organization}/${repository}/pullrequests/${pr}/approve")

    def HttpURLConnection conn = (HttpURLConnection)url.openConnection()
    conn.setRequestProperty( "Content-Type", "application/x-www-form-urlencoded")
    conn.setRequestProperty( "Authorization", config.basicAuth)
    conn.setRequestMethod("DELETE")

    try {
        statusCode = conn.getResponseCode()
        if (statusCode==204) {
            println "[Unapprove RP] PR succesfully unapproved"
            prUnapproved = true
        } else if (statusCode==404) {
            println "[Unapprove RP] PR not approved by this user"
            prUnapproved = true
        } else {
            println "[Unapprove PR] Connection status code: ${statusCode}"
            println "[Unapprove PR] URL ${url.toString()}"
            println "[Unapprove PR] Server response:"
            response = displayServerResponse(conn.getErrorStream())
            println "[Unapprove PR] ${response}"
            println "[Unapprove PR] -----"
        }
    } catch (Exception e) {
        println "[Unapprove PR] Error connecting to ${url.toString()}"
        println e.getMessage()
    } finally {
        if (conn != null) {
            conn.disconnect();
        }
    }
    return prUnapproved
}
